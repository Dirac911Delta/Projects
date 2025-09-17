# Laser Dosage Time Calculator

## Purpose

This Excel workbook calculates the required laser exposure time based on a desired dose (energy per unit area) and detector specifications. It's designed for applications where you need to determine how long to expose a detector to laser light to achieve a specific energy dose.

## How to Use the Calculator

### 1. Setting Input Parameters

Open the **Inputs** sheet and modify the following values:

- **Desired Dose (B1)**: Enter the target energy dose in J/cm² (default: 5 J/cm²)
- **Detector Diameter (B2)**: Enter the detector diameter in millimeters (default: 19 mm)

The following values are automatically calculated:
- **Detector Area (B3)**: Automatically computed from the diameter
- **Constant factor k (B8)**: Pre-calculated factor for efficiency

### 2. Viewing Results

Switch to the **Data** sheet to see:

- **Measurement Table**: Shows power levels, currents, and calculated exposure times
- **Chart**: Visual representation of the relationship between laser power and required exposure time

### 3. Interpreting Results

For each power level in the table:
- **Required Time (s)**: Time needed in seconds to achieve the desired dose
- **Required Time (min)**: Same time converted to minutes for convenience

## Formula and Unit Assumptions

The calculator uses the fundamental relationship:

```
t[s] = 1000 * D(J/cm²) * A(cm²) / P(mW)
```

Where:
- **t** = Required exposure time in seconds
- **D** = Desired dose in J/cm²
- **A** = Detector area in cm²
- **P** = Laser power in milliwatts

### Unit Conversions
- Detector diameter is entered in **millimeters** but converted to cm for area calculation
- Area formula: A = π × (diameter_cm/2)² = π × (diameter_mm/20)²
- The factor of 1000 converts between J (joules) and mW·s (milliwatt-seconds)

## Chart Updates

The scatter plot chart automatically updates when you change the input parameters (desired dose or detector diameter). This allows you to visualize how different settings affect the required exposure times across various power levels.

## Sample Verification Results

With the default settings (Dose = 5 J/cm², Diameter = 19 mm):
- Detector area = 2.835 cm²
- Expected results:
  - 33 mW → 429.59 s (7.16 min)
  - 50 mW → 283.53 s (4.73 min)
  - 70 mW → 202.52 s (3.38 min)
  - 90 mW → 157.52 s (2.63 min)
  - 110 mW → 128.88 s (2.15 min)

## Data Validation

The workbook includes input validation:
- **Desired Dose**: Must be a positive number
- **Detector Diameter**: Must be a positive number

## Sheet Protection

The **Inputs** sheet is protected to prevent accidental modification of formulas, but the input cells (B1 and B2) remain editable. If you need to make other changes, you can unprotect the sheet (password: "laser123").

## Technical Notes

- The workbook uses Excel tables for structured data management
- Formulas use absolute references to ensure proper calculation when inputs change
- Number formatting is applied for readability (3 decimal places for area, 2 for times)
- The chart uses dynamic ranges that automatically include all measurement data

## Troubleshooting

If calculations don't update automatically:
1. Check that Excel calculation is set to "Automatic" (Formulas > Calculation Options)
2. Press F9 to force recalculation
3. Ensure the input values are valid numbers greater than zero

For technical support or modifications, refer to the source Python script `create_laser_workbook.py` included in this repository.