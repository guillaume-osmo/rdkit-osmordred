#!/usr/bin/env python3
"""
Script to apply the Gasteiger NaN patch to calcRNCG_RPCG function.
"""

import re

# Read the file
with open('Code/GraphMol/Descriptors/Osmordred.cpp', 'r') as f:
    content = f.read()

# Find the function using a more flexible pattern
# Look for the function signature and the computeGasteigerCharges call
pattern = r'(std::vector<double> calcRNCG_RPCG\(const RDKit::ROMol& mol\)\{[^}]*?computeGasteigerCharges\(\*hmol, 12, true\);[^}]*?return \{maxneg/totalneg, maxpos/totalpos\};\s*\n\s*\n\})'

# Try a simpler approach: find the function start and end
func_start = content.find('std::vector<double> calcRNCG_RPCG(const RDKit::ROMol& mol){')
if func_start == -1:
    print("ERROR: Function not found!")
    exit(1)

# Find the matching closing brace
brace_count = 0
func_end = func_start
for i in range(func_start, len(content)):
    if content[i] == '{':
        brace_count += 1
    elif content[i] == '}':
        brace_count -= 1
        if brace_count == 0:
            func_end = i + 1
            break

if func_end == func_start:
    print("ERROR: Could not find function end!")
    exit(1)

# Extract the original function
original_func = content[func_start:func_end]

# Create the patched version
patched_func = '''std::vector<double> calcRNCG_RPCG(const RDKit::ROMol& mol){
        // subclass of CPSA using only 2D descriptors available in Mordred v1
         std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));

        double maxpos = 0;
        double maxneg = 0;
        double totalneg =0;
        double totalpos =0;
        bool gasteiger_failed = false;

        // Try to compute Gasteiger charges, catch any exceptions
        try {
            computeGasteigerCharges(*hmol, 12, true);
        } catch (...) {
            // Gasteiger failed (e.g., Hg not in parameters) - mark as failed
            gasteiger_failed = true;
            // Set all charges to NaN
            for (auto &atom : hmol->atoms()) {
                atom->setProp(common_properties::_GasteigerCharge, std::numeric_limits<double>::quiet_NaN());
            }
        }

        // If Gasteiger failed, return NaN for both values
        if (gasteiger_failed) {
            return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        }

        // Check if any charges are NaN (even if computeGasteigerCharges didn't throw)
        bool has_nan_charges = false;
        for (auto &atom : hmol->atoms()) {
            double ch = atom->getProp<double>(common_properties::_GasteigerCharge);

            if (std::isnan(ch)) {
                has_nan_charges = true;
                break;
            }

            if (atom->hasProp(common_properties::_GasteigerHCharge)) {
                double hch = atom->getProp<double>(common_properties::_GasteigerHCharge);
                if (std::isnan(hch)) {
                    has_nan_charges = true;
                    break;
                }
                ch += hch;
            }

            // Check again after adding H charge
            if (std::isnan(ch)) {
                has_nan_charges = true;
                break;
            }

            if (ch < 0) {
                totalneg += -ch;
                if (-ch > maxneg) {
                    maxneg = -ch;
                }
            }
            else if ( ch > 0) {
                totalpos += ch;
                if (ch > maxpos) {
                    maxpos = ch;
                }
            }
        }

        // If any charges were NaN, return NaN for both values
        if (has_nan_charges) {
            return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        }

        if (totalneg == 0 || totalpos == 0){
            return {0.,0.};
        }
        return {maxneg/totalneg, maxpos/totalpos};


}'''

# Replace in content
new_content = content[:func_start] + patched_func + content[func_end:]

# Write back
with open('Code/GraphMol/Descriptors/Osmordred.cpp', 'w') as f:
    f.write(new_content)

print("✅ Patch applied successfully!")
print(f"   Function replaced: {len(original_func)} chars -> {len(patched_func)} chars")

