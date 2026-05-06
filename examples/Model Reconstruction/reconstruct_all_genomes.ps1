[CmdletBinding()]
param(
    [string]$Python = "python",
    [string]$Uv = "uv",
    [string]$GenomeDir,
    [string]$OutputDir,
    [string]$DatabasePath,
    [string]$TemplateDir,
    [string]$FastaGlob,
    [ValidateSet("rast", "fasta-description", "none")]
    [string]$AnnotationSource = "rast",
    [switch]$DryRunTemplateSelection,
    [switch]$ForceTemplateBuild,
    [switch]$InstallMissingPackages,
    [switch]$AddAtpm,
    [switch]$AllowAllNonGprReactions,
    [switch]$ClassicBiomass,
    [double]$Gc = 0.5
)

$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $PSCommandPath
$RepoRoot = (Resolve-Path (Join-Path $ScriptDir "..\..")).Path

function Get-FullPath {
    param([Parameter(Mandatory = $true)][string]$PathValue)

    if ([System.IO.Path]::IsPathRooted($PathValue)) {
        return [System.IO.Path]::GetFullPath($PathValue)
    }
    return [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $PathValue))
}

if (-not $GenomeDir) {
    $GenomeDir = Join-Path $RepoRoot "genomes"
} else {
    $GenomeDir = Get-FullPath $GenomeDir
}

if (-not $OutputDir) {
    $OutputDir = Join-Path $RepoRoot "reconstructed_models"
} else {
    $OutputDir = Get-FullPath $OutputDir
}

if (-not $DatabasePath) {
    $DatabasePath = Join-Path $RepoRoot "ModelSEEDDatabase"
} else {
    $DatabasePath = Get-FullPath $DatabasePath
}

if (-not $TemplateDir) {
    $TemplateDir = Join-Path $OutputDir "templates"
} else {
    $TemplateDir = Get-FullPath $TemplateDir
}

$BuildTemplateScript = Join-Path $ScriptDir "build_template_from_modelseed_database.py"
$ReconstructScript = Join-Path $ScriptDir "reconstruct_non_gapfilled_model.py"

if (-not (Test-Path -LiteralPath $GenomeDir -PathType Container)) {
    throw "Genome directory not found: $GenomeDir"
}
if (-not (Test-Path -LiteralPath $DatabasePath -PathType Container)) {
    throw "ModelSEEDDatabase directory not found: $DatabasePath"
}
if (-not (Test-Path -LiteralPath $BuildTemplateScript -PathType Leaf)) {
    throw "Template builder script not found: $BuildTemplateScript"
}
if (-not (Test-Path -LiteralPath $ReconstructScript -PathType Leaf)) {
    throw "Reconstruction script not found: $ReconstructScript"
}

$NucleotideInputs = @(
    Get-ChildItem -LiteralPath $GenomeDir -File -Filter "*.fna" -ErrorAction SilentlyContinue
    Get-ChildItem -LiteralPath $GenomeDir -File -Filter "*.fna.gz" -ErrorAction SilentlyContinue
)

if ($NucleotideInputs.Count -gt 0 -and -not $DryRunTemplateSelection) {
    & $Python -c "import pyrodigal" 2>$null
    if ($LASTEXITCODE -ne 0) {
        if ($InstallMissingPackages) {
            Write-Host "Installing pyrodigal for nucleotide FASTA translation"
            & $Uv pip install pyrodigal
            if ($LASTEXITCODE -ne 0) {
                throw "Failed to install pyrodigal"
            }
        } else {
            throw "Found .fna/.fna.gz inputs, but pyrodigal is not installed. Run `uv pip install pyrodigal` or rerun this script with -InstallMissingPackages."
        }
    }
}

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
New-Item -ItemType Directory -Force -Path $TemplateDir | Out-Null

$TemplateSpecs = @(
    @{
        Template = "GramPositive"
        Id       = "GramPositive.latest.modeltemplate"
        File     = "GramPositive.latest.modeltemplate.json"
    },
    @{
        Template = "GramNegative"
        Id       = "GramNegative.latest.modeltemplate"
        File     = "GramNegative.latest.modeltemplate.json"
    }
)

foreach ($Spec in $TemplateSpecs) {
    $Target = Join-Path $TemplateDir $Spec.File
    if ($ForceTemplateBuild -or -not (Test-Path -LiteralPath $Target -PathType Leaf)) {
        Write-Host "Building updated $($Spec.Template) template -> $Target"
        & $Python $BuildTemplateScript `
            --database-path $DatabasePath `
            --template $Spec.Template `
            --template-id $Spec.Id `
            --name $Spec.Id `
            --biochem-ref "local/ModelSEEDDatabase" `
            --output $Target
        if ($LASTEXITCODE -ne 0) {
            throw "Template build failed for $($Spec.Template)"
        }
    } else {
        Write-Host "Using existing updated $($Spec.Template) template: $Target"
    }
}

$ReconstructArgs = @(
    $ReconstructScript,
    "--fasta-dir", $GenomeDir,
    "--output-dir", $OutputDir,
    "--updated-template-dir", $TemplateDir,
    "--biochem-path", $DatabasePath,
    "--annotation-source", $AnnotationSource,
    "--gc", $Gc.ToString([System.Globalization.CultureInfo]::InvariantCulture)
)

if ($FastaGlob) {
    $ReconstructArgs += @("--fasta-glob", $FastaGlob)
}
if ($DryRunTemplateSelection) {
    $ReconstructArgs += "--dry-run-template-selection"
}
if ($AddAtpm) {
    $ReconstructArgs += "--add-atpm"
}
if ($AllowAllNonGprReactions) {
    $ReconstructArgs += "--allow-all-non-gpr-reactions"
}
if ($ClassicBiomass) {
    $ReconstructArgs += "--classic-biomass"
}

if ($FastaGlob) {
    Write-Host "Running reconstruction over $GenomeDir ($FastaGlob)"
} else {
    Write-Host "Running reconstruction over $GenomeDir (*.faa, *.faa.gz, *.fna, *.fna.gz)"
}
& $Python @ReconstructArgs
if ($LASTEXITCODE -ne 0) {
    throw "Batch reconstruction failed"
}
