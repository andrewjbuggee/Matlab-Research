# Quality Control Decision Tree

## Visual Guide: How the Function Processes Each Sounding

```
┌─────────────────────────────────────┐
│  Start Processing Sounding          │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│  Parse Header Information            │
│  • Data type, site, location        │
│  • Release time, quality info       │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│  Read All Data Rows                  │
│  • Variable names, units            │
│  • Data matrix (N × 21)             │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│  Check Each Row for Errors          │
│                                     │
│  Erroneous if ANY of:               │
│  • Pressure ≥ 9999 mb              │
│  • Temperature ≥ 999°C             │
│  • Temperature ≤ -999°C            │
│  • Dewpoint ≥ 999°C                │
│  • Dewpoint ≤ -999°C               │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│  Count Erroneous Rows               │
│  Calculate: % = (bad/total) × 100  │
└──────────────┬──────────────────────┘
               │
               ▼
        ┌──────┴──────┐
        │             │
    % > 50%       % ≤ 50%
        │             │
        ▼             ▼
┌───────────────┐  ┌───────────────────┐
│ SKIP SOUNDING │  │ CLEAN SOUNDING    │
│               │  │                   │
│ • Print       │  │ • Remove bad rows │
│   warning     │  │ • Print removal   │
│ • Return      │  │   count & %       │
│   empty []    │  │ • Return cleaned  │
│ • Don't add   │  │   sounding        │
│   to output   │  │ • Add to output   │
└───────────────┘  └───────┬───────────┘
        │                  │
        └────────┬─────────┘
                 │
                 ▼
        ┌────────────────┐
        │  Next Sounding │
        └────────────────┘
```

## Example Scenarios

### Scenario 1: Good Sounding (5% bad)
```
Input:  100 rows total
        5 rows with pressure = 9999
        
Check:  5/100 = 5% erroneous

Decision: 5% ≤ 50% ✓
          
Action: Remove 5 rows
        Keep sounding with 95 rows
        
Output: ✓ Included in structure
```

### Scenario 2: Marginal Sounding (48% bad)
```
Input:  90 rows total
        44 rows with dewpoint = 999
        
Check:  44/90 = 48.9% erroneous

Decision: 48.9% ≤ 50% ✓
          
Action: Remove 44 rows
        Keep sounding with 46 rows
        
Output: ✓ Included in structure
```

### Scenario 3: Poor Sounding (68% bad)
```
Input:  100 rows total
        68 rows with temperature = 999
        
Check:  68/100 = 68% erroneous

Decision: 68% > 50% ✗
          
Action: Skip entire sounding
        Warning message printed
        
Output: ✗ NOT included in structure
```

## Output Structure Impact

### File Contains: 15 Soundings

```
Sounding 1:  5% bad  → ✓ Included (index 1)
Sounding 2:  3% bad  → ✓ Included (index 2)
Sounding 3: 48% bad  → ✓ Included (index 3)
Sounding 4:  0% bad  → ✓ Included (index 4)
Sounding 5: 12% bad  → ✓ Included (index 5)
Sounding 6: 68% bad  → ✗ SKIPPED
Sounding 7:  1% bad  → ✓ Included (index 6)
Sounding 8: 15% bad  → ✓ Included (index 7)
Sounding 9:  8% bad  → ✓ Included (index 8)
Sounding 10: 22% bad → ✓ Included (index 9)
Sounding 11: 75% bad → ✗ SKIPPED
Sounding 12:  0% bad → ✓ Included (index 10)
Sounding 13:  4% bad → ✓ Included (index 11)
Sounding 14: 90% bad → ✗ SKIPPED
Sounding 15:  6% bad → ✓ Included (index 12)

Result: data array has 12 elements (not 15!)
```

### Accessing Results:
```matlab
% DON'T assume original sounding numbers
for original_num = 1:15  % ✗ WRONG - may not exist
    process(data(original_num));
end

% DO use actual structure length
for i = 1:length(data)  % ✓ CORRECT
    process(data(i));
end
```

## Console Output Example

```
Found 15 sounding(s) in file: VOCALS_2008_5mb_20081025.cls

--- Processing Sounding 1/15 ---
  Removed 5 erroneous data row(s) from sounding 1 (5.0% of data)
Release Site: SCQN/85575
Number of data points: 95

--- Processing Sounding 2/15 ---
Release Site: R/V Ron Brown
Number of data points: 200

--- Processing Sounding 3/15 ---
  Removed 44 erroneous data row(s) from sounding 3 (48.9% of data)
Release Site: SCQN/85575
Number of data points: 46

...

--- Processing Sounding 6/15 ---
  WARNING: 68.0% of data is erroneous (>50%) - SKIPPING this sounding
  Sounding 6 skipped due to excessive erroneous data

--- Processing Sounding 7/15 ---
Release Site: R/V Jose Olaya
Number of data points: 195

...

=== Successfully read 12 out of 15 sounding(s) ===
```

## Quality Assurance Flowchart

```
┌──────────────────────────────────────────┐
│ User calls: data = read_cls_file(file)  │
└─────────────────┬────────────────────────┘
                  │
                  ▼
         ┌────────────────┐
         │ Find soundings │
         └────────┬───────┘
                  │
          ┌───────┴────────┐
          │ For each       │
          │ sounding:      │
          └───────┬────────┘
                  │
    ┌─────────────┼─────────────┐
    │             │             │
    ▼             ▼             ▼
┌────────┐  ┌─────────┐  ┌─────────┐
│ Parse  │  │ Check   │  │ Apply   │
│ header │→ │ quality │→ │ filter  │
└────────┘  └─────────┘  └────┬────┘
                               │
                      ┌────────┴────────┐
                      │                 │
                  Good (≤50%)      Poor (>50%)
                      │                 │
                      ▼                 ▼
                ┌──────────┐      ┌─────────┐
                │ Add to   │      │ Skip    │
                │ output   │      │ (empty) │
                └────┬─────┘      └────┬────┘
                     │                 │
                     └────────┬────────┘
                              │
                              ▼
                    ┌─────────────────┐
                    │ Next sounding   │
                    └─────────┬───────┘
                              │
                              ▼
                    ┌─────────────────┐
                    │ Return data     │
                    │ (validated)     │
                    └─────────────────┘
```

## Summary Statistics

After processing, you get:

```matlab
% Total soundings found in file
fprintf('File contained: %d soundings\n', num_found);

% How many passed quality control
fprintf('Quality-validated: %d soundings\n', length(data));

% How many were excluded
fprintf('Excluded: %d soundings (>50%% bad data)\n', 
        num_found - length(data));
```

## Key Takeaways

1. 🔍 **Every row checked** for pressure, temperature, AND dewpoint errors
2. 📊 **Percentage calculated** to determine sounding quality
3. ⚖️ **50% threshold** separates good from poor soundings
4. ✂️ **Bad rows removed** from good soundings
5. 🗑️ **Poor soundings excluded** entirely
6. ✅ **Output guaranteed** to be high quality
7. 📢 **Clear feedback** about what was done

**Result: Only reliable atmospheric profiles in your dataset!**
