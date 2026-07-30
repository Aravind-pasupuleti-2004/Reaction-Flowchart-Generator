# AI Chemical Structure Generator

A completely local desktop application that automatically generates accurate 2D molecular structures from SMILES, IUPAC names, InChI, or common compound names. No cloud services or paid APIs required.

## Features

- **Single Compound Mode**: Enter one compound at a time via SMILES, IUPAC name, common name, or InChI
- **Batch Processing**: Upload Excel/CSV files with multiple compounds
- **High-quality Rendering**: PNG and SVG output
- **Compound Information**: Molecular formula, weight, SMILES, InChI, atom/ring counts
- **Download Options**: PNG, SVG, Molfile, batch report (Excel), images ZIP

## Installation

1. Install Python 3.12+

2. Install dependencies:

```bash
cd chemical_structure_generator
pip install -r requirements.txt
```

> **Note**: If the above doesn't install RDKit, try:
> ```bash
> pip install rdkit-pypi
> ```

3. Install Java Runtime (JRE) — required for IUPAC name conversion via pyopsin:
   Download from https://adoptium.net/ or run:
   ```bash
   winget install EclipseAdoptium.Temurin.21.JRE
   ```
   Set `JAVA_HOME` environment variable to the JRE install path.

4. Run the application:

```bash
streamlit run app.py
```

## Usage

### Single Compound
1. Select "Single Compound" mode
2. Enter SMILES, IUPAC name, common name, or InChI
3. Click "Generate"
4. View structure and download as PNG/SVG/Molfile

### Batch Processing
1. Select "Batch Processing" mode
2. Upload an Excel (.xlsx) or CSV file
3. Click "Generate All Structures"
4. Download the report and images ZIP

## Input Types

| Type | Example |
|------|---------|
| SMILES | `CCO` |
| IUPAC | `ethanol` |
| Common name | `aspirin` |
| InChI | `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3` |

## Project Structure

```
chemical_structure_generator/
├── app.py              # Streamlit UI
├── converter.py        # Input detection & conversion logic
├── renderer.py         # RDKit molecule rendering (PNG, SVG)
├── batch_processor.py  # Excel/CSV batch processing
├── pubchem_lookup.py   # Optional PubChem API lookup
├── utils.py            # Utility functions
├── outputs/            # Generated outputs
│   ├── images/
│   ├── reports/
│   └── zip/
├── requirements.txt
└── README.md
```

## Requirements

- Python 3.12+
- Java Runtime (JRE) 17+
- RDKit
- Streamlit
- Pandas, OpenPyXL
- pyopsin (for offline IUPAC name conversion via OPSIN)
- Pillow

## License

MIT
