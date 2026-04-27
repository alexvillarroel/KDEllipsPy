#!/bin/bash

# Example script: How to set up and run an event-run folder

# Step 1: Create run folder structure
echo "Creating example run folder structure..."
RUN_DIR="example_event_run"
mkdir -p "$RUN_DIR/DATA"

# Step 2: Copy configuration (customize as needed)
echo "Copying example configuration..."
cp input.ctl "$RUN_DIR/"

# Step 3: Copy or link observed data
echo "Linking observed data..."
ln -sf ../DATA/real_disp_x "$RUN_DIR/DATA/"
ln -sf ../DATA/real_disp_y "$RUN_DIR/DATA/"
ln -sf ../DATA/real_disp_z "$RUN_DIR/DATA/"

# Step 4: Run inversion
echo ""
echo "Running inversion..."
python main.py --run "$RUN_DIR"

# Step 5: Check results
echo ""
echo "Results saved in: $RUN_DIR/output/"
echo ""
ls -lh "$RUN_DIR/output/"

# Step 6: View results
echo ""
echo "Best model summary:"
cat "$RUN_DIR/output/na_results.json" | jq '.best_model'
