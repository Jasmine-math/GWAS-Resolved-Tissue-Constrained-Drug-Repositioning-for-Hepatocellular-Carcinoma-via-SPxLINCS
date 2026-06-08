# ====== Configure paths ======
$Folder  = "D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\top_drug keytissue"
$NodeOut = "D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\cyto_nodes_drug_gwas_final.csv"
$EdgeOut = "D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\cyto_edges_drug_gwas.csv"

# Ensure output directory exists
$OutDir = Split-Path -Path $NodeOut -Parent
if (-not (Test-Path -LiteralPath $OutDir)) {
    New-Item -ItemType Directory -Path $OutDir -Force | Out-Null
}

# ====== Plain CSV writer (UTF-8 without BOM, minimal quoting) ======
function Export-PlainCsv {
    param(
        [Parameter(Mandatory = $true)]
        [System.Collections.IEnumerable]$Rows,

        [Parameter(Mandatory = $true)]
        [string]$Path
    )

    $rowList = @($Rows)

    if ($rowList.Count -eq 0) {
        [System.IO.File]::WriteAllText($Path, "", [System.Text.UTF8Encoding]::new($false))
        return
    }

    $headers = @($rowList[0].PSObject.Properties.Name)
    $sb = New-Object System.Text.StringBuilder

    # header
    [void]$sb.AppendLine(($headers -join ","))

    foreach ($row in $rowList) {
        $fields = foreach ($h in $headers) {
            $v = $row.$h

            if ($null -eq $v) {
                ""
            } else {
                $s = $v.ToString()
                $s = $s -replace "`r", " " -replace "`n", " "
                $s = $s.Trim()

                if ($s.Contains(",") -or $s.Contains('"')) {
                    $s = $s -replace '"', '""'
                    '"' + $s + '"'
                } else {
                    $s
                }
            }
        }

        [void]$sb.AppendLine(($fields -join ","))
    }

    [System.IO.File]::WriteAllText($Path, $sb.ToString(), [System.Text.UTF8Encoding]::new($false))
}

# ====== Auto-detect delimiter (CSV or TSV) and import ======
function Import-WithAutoDelimiter {
    param([string]$Path)

    $firstLine  = Get-Content -LiteralPath $Path -TotalCount 1
    $tabCount   = ([regex]::Matches($firstLine, "`t")).Count
    $commaCount = ([regex]::Matches($firstLine, ",")).Count

    if ($tabCount -gt $commaCount) {
        $delimiter = "`t"
    } else {
        $delimiter = ","
    }

    Import-Csv -LiteralPath $Path -Delimiter $delimiter
}

# ====== Initialize containers ======
$files       = Get-ChildItem -Path $Folder -Filter *.csv -File
$index       = @{}   # $index[drugLower][gwasId] = @{ minq = double; count = int }
$drugNameMap = @{}   # preserve first-seen original drug casing
$gwasSet     = @{}   # set of GWAS IDs

$culture = [System.Globalization.CultureInfo]::InvariantCulture

# ====== Iterate input files ======
foreach ($f in $files) {

    # Parse GWAS ID from filename
    # substring between 'fallback3_' and '_hg38'; fallback to basename
    $name = $f.BaseName
    $m = [regex]::Match(
        $name,
        'fallback3_(.+)_hg38$',
        [System.Text.RegularExpressions.RegexOptions]::IgnoreCase
    )

    if ($m.Success) {
        $gwasId = $m.Groups[1].Value.Trim()
    } else {
        $gwasId = $name.Trim()
    }

    $gwasId = ($gwasId -replace '\s+', '_').Trim()
    if ($gwasId -eq "") { continue }

    if (-not $gwasSet.ContainsKey($gwasId)) {
        $gwasSet[$gwasId] = $true
    }

    $df = Import-WithAutoDelimiter -Path $f.FullName
    if ($null -eq $df -or $df.Count -eq 0) {
        continue
    }

    $cols = $df[0].PSObject.Properties.Name

    # Find Drug column (prefer exact 'Drug'; otherwise use first column)
    $drugCol = $null
    foreach ($c in $cols) {
        if ($c -match '^(?i)Drug$') {
            $drugCol = $c
            break
        }
    }
    if (-not $drugCol) {
        $drugCol = $cols[0]
    }

    # Find q column
    $qCol = $null
    foreach ($c in $cols) {
        if ($c -match '^(?i)q[\s_-]*value$' -or $c -match '^(?i)q$') {
            $qCol = $c
            break
        }
    }

    if (-not $qCol) {
        foreach ($c in $cols) {
            if ( ($c -replace '[\s_-]', '') -match '^(?i)qvalue$' ) {
                $qCol = $c
                break
            }
        }
    }

    if (-not $qCol) {
        Write-Warning "Skip (q column not found): $($f.Name)"
        continue
    }

    foreach ($r in $df) {
        $drug = $r.$drugCol
        if ($null -eq $drug) { continue }

        $drug = $drug.ToString().Trim()
        if ($drug -eq "") { continue }

        $qStr = $r.$qCol
        if ($null -eq $qStr) { continue }

        $qStr = $qStr.ToString().Trim()
        if ($qStr -eq "") { continue }

        $qVal = 0.0
        if (-not [double]::TryParse(
            $qStr,
            [System.Globalization.NumberStyles]::Float,
            $culture,
            [ref]$qVal
        )) {
            continue
        }

        if ($qVal -lt 0) { continue }
        if ($qVal -lt [double]::Epsilon) {
            $qVal = [double]::Epsilon
        }

        $key = $drug.ToLowerInvariant()

        if (-not $index.ContainsKey($key)) {
            $index[$key] = @{}
            $drugNameMap[$key] = $drug
        }

        if (-not $index[$key].ContainsKey($gwasId)) {
            $index[$key][$gwasId] = @{
                minq  = $qVal
                count = 1
            }
        } else {
            if ($qVal -lt $index[$key][$gwasId].minq) {
                $index[$key][$gwasId].minq = $qVal
            }
            $index[$key][$gwasId].count = $index[$key][$gwasId].count + 1
        }
    }
}

# ====== Prepare sorted lists ======
$gwasList = $gwasSet.Keys | Sort-Object
$drugKeys = $index.Keys | Sort-Object

# ====== Build final node table ======
$nodeRows = New-Object System.Collections.Generic.List[object]

# Drug nodes
foreach ($key in $drugKeys) {
    $drugName = $drugNameMap[$key].Trim()

    $props = New-Object 'System.Collections.Specialized.OrderedDictionary'
    $props['name']            = $drugName
    $props['shared_name']     = $drugName
    $props['node_type']       = 'drug'
    $props['label_text']      = $drugName
    $props['traitname']       = $drugName

    $pieTotal = 0.0
    foreach ($g in $gwasList) {
        if ($index[$key].ContainsKey($g)) {
            $minq = $index[$key][$g].minq
            $neglog = -([math]::Log10($minq))
            $neglogRounded = [math]::Round($neglog, 6)
            $props[$g] = $neglogRounded
            $pieTotal += $neglogRounded
        } else {
            $props[$g] = $null
        }
    }

    $props['pie_total_neglog'] = [math]::Round($pieTotal, 6)
    $nodeRows.Add((New-Object psobject -Property $props))
}

# GWAS nodes
foreach ($g in $gwasList) {
    $g2 = $g.Trim()

    $props = New-Object 'System.Collections.Specialized.OrderedDictionary'
    $props['name']            = $g2
    $props['shared_name']     = $g2
    $props['node_type']       = 'gwas'
    $props['label_text']      = $g2
    $props['traitname']       = $g2

    foreach ($gg in $gwasList) {
        $props[$gg] = $null
    }

    $props['pie_total_neglog'] = $null
    $nodeRows.Add((New-Object psobject -Property $props))
}

# ====== Build edge table ======
$edgeRows = New-Object System.Collections.Generic.List[object]

foreach ($key in $drugKeys) {
    $drugName = $drugNameMap[$key].Trim()

    foreach ($g in ($index[$key].Keys | Sort-Object)) {
        $g2 = $g.Trim()
        $minq = $index[$key][$g].minq
        $neglog = -([math]::Log10($minq))
        $count = $index[$key][$g].count

        $edgeRows.Add([pscustomobject]@{
            source         = $drugName
            target         = $g2
            interaction    = 'Drug-GWAS'
            count          = $count
            min_q          = $minq
            neglog10_min_q = [math]::Round($neglog, 6)
        })
    }
}

# ====== Export using plain CSV ======
Export-PlainCsv -Rows $nodeRows -Path $NodeOut
Export-PlainCsv -Rows $edgeRows -Path $EdgeOut

# ====== Summary ======
Write-Host "Generated final node table:" $NodeOut
Write-Host "Generated edge table:" $EdgeOut
Write-Host "Drugs:" $drugKeys.Count " GWAS:" $gwasList.Count " Edges:" $edgeRows.Count