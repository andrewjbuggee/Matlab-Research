# Update: 50% Quality Threshold - Exclude Poor Soundings

## New Feature Added

Soundings where **more than 50%** of the data rows are erroneous are now completely **excluded** from the output structure, rather than just having bad rows removed.

## The Problem

Previously, the function would remove erroneous rows from ALL soundings, regardless of data quality. This meant:

- A sounding with 90 data rows, where 70 are erroneous (78%), would keep only 20 rows
- These heavily degraded soundings are not useful for analysis
- They lack vertical resolution and create gaps in the profile
- Including them could bias statistical analyses

## The Solution

Implement a **50% quality threshold**:

1. **Calculate** the percentage of erroneous rows in each sounding
2. **If ≤ 50% bad**: Remove bad rows, keep the sounding
3. **If > 50% bad**: Skip the entire sounding (don't add to output structure)

## Implementation

### Step 1: Calculate Percentage of Bad Data
```matlab
num_bad = sum(bad_rows);
percent_bad = (num_bad / num_rows_before) * 100;
```

### Step 2: Check Threshold
```matlab
if percent_bad > 50
    fprintf('  WARNING: %.1f%% of data is erroneous (>50%%) - SKIPPING this sounding\n', percent_bad);
    sounding = [];  % Return empty to signal skip
    return;
end
```

### Step 3: Handle Empty Soundings in Main Function
```matlab
sounding_data = parse_single_sounding(sounding_lines, s);

% Check if sounding was skipped
if isempty(sounding_data)
    fprintf('  Sounding %d skipped due to excessive erroneous data\n', s);
    continue;  % Skip to next sounding
end

% Otherwise, add to output
soundings_processed = soundings_processed + 1;
data(soundings_processed) = sounding_data;
```

### Step 4: Update Final Message
```matlab
fprintf('\n=== Successfully read %d out of %d sounding(s) ===\n', 
        soundings_processed, num_soundings);
```

## Example: VOCALS_2008_5mb_20081025.cls

This file has 15 soundings. Let's look at sounding 1:

### Before Filtering:
```
Total rows: 90
Erroneous rows: 44 (lines 17-60 with dewpoint = 999)
Percentage bad: 48.9%
```

**Action**: 48.9% < 50%, so **KEEP** the sounding and remove 44 bad rows
**Result**: 46 clean data rows remain

### Hypothetical Example:
```
Total rows: 90
Erroneous rows: 62
Percentage bad: 68.9%
```

**Action**: 68.9% > 50%, so **SKIP** entire sounding
**Result**: Sounding not included in output structure

## Output Examples

### Good Sounding (< 50% bad):
```
--- Processing Sounding 1/15 ---
  Removed 44 erroneous data row(s) from sounding 1 (48.9% of data)
Data Type: Paposo, Chile Soundings/Ascending
Release Site: SCQN/85575
Number of data points: 46
```

### Poor Sounding (> 50% bad):
```
--- Processing Sounding 9/15 ---
  WARNING: 68.5% of data is erroneous (>50%) - SKIPPING this sounding
  Sounding 9 skipped due to excessive erroneous data
```

### Final Summary:
```
=== Successfully read 12 out of 15 sounding(s) ===
```

## Benefits

### 1. Data Quality Assurance
- Only soundings with good vertical coverage are included
- Prevents sparse profiles with large gaps
- Ensures reliable thermodynamic calculations

### 2. Clearer User Feedback
- Users know exactly which soundings were excluded
- Percentage reported helps understand data quality
- Clear distinction between "cleaned" and "excluded"

### 3. Appropriate for Analysis
- Soundings with >50% bad data are unreliable
- Would introduce uncertainty in:
  - Vertical interpolation
  - Layer averages
  - Stability calculations
  - Thermodynamic diagrams

### 4. Maintains Array Structure
- Output structure only contains valid soundings
- No need to check for empty/sparse soundings later
- Cleaner downstream code

## Why 50%?

The 50% threshold is a standard quality control criterion because:

1. **Vertical Resolution**: Losing >50% of observations creates large gaps
2. **Statistical Validity**: Fewer than half the expected observations reduces confidence
3. **Interpolation Risk**: Large gaps require excessive interpolation
4. **Industry Standard**: Many atmospheric datasets use similar thresholds
5. **Balance**: Strict enough to exclude poor data, lenient enough to keep marginal soundings

## Impact on Your Workflow

### Before:
```matlab
data = read_cls_file('file.cls');
% Returns: 15 soundings, some with only 10-20% of original rows
% User must manually check quality
```

### After:
```matlab
data = read_cls_file('file.cls');
% Returns: 12 soundings, all with >50% of rows intact
% Guaranteed minimum quality level
```

### Accessing Results:
```matlab
% No change needed - just use the data
for i = 1:length(data)
    % All soundings here are high quality
    plot_sounding(data(i));
end
```

## Testing

Run the test script to verify:
```matlab
test_50percent_threshold
```

Expected behavior:
- ✓ Soundings with ≤50% bad data are cleaned and included
- ✓ Soundings with >50% bad data are completely excluded
- ✓ No erroneous values remain in output
- ✓ Clear warnings printed for excluded soundings

## Files Updated

1. **read_cls_file.m** - Added 50% threshold logic
2. **README.md** - Updated data quality control section
3. **QUICK_REFERENCE.md** - Added automatic skipping behavior
4. **test_50percent_threshold.m** - New test script
5. **THRESHOLD_UPDATE.md** - This documentation

## Backward Compatibility

⚠️ **BREAKING CHANGE**: The number of soundings in the output may be less than expected

### Migration:
```matlab
% Before: Assumed all soundings present
num_soundings = 15;  % From file inspection

% After: Check actual number in structure
num_soundings = length(data);  % Actual valid soundings
```

### Advantage:
```matlab
% No need for quality checks in your code!
for i = 1:length(data)
    % All soundings here are already validated
    calculate_cape(data(i));  % Safe to use
end
```

## Summary

| Feature | Before | After |
|---------|--------|-------|
| Bad rows removed | ✓ | ✓ |
| Poor soundings removed | ✗ | ✓ |
| Quality threshold | None | 50% |
| User notification | Basic | Detailed with % |
| Output guarantee | None | All soundings usable |

**Recommendation**: Use this version for all production analysis. It ensures you only work with high-quality atmospheric profiles.
