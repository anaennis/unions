#!/usr/bin/env python3
"""
Update dolphot.param file with values from CONFIG.py
"""

import os
import sys

# Import CONFIG
import CONFIG

def update_dolphot_param(filename):
    """Update dolphot.param file with values from CONFIG.py"""
    
    if not os.path.exists(filename):
        print(f"ERROR: {filename} not found!")
        return False
    
    # Mapping from CONFIG to dolphot.param keys
    updates = {
        'Nimg': str(CONFIG.Nimg),
        'img_RAper': str(CONFIG.img_RAPER),
        'img_RPSF': str(CONFIG.img_RPSF),
        'img_aprad': str(CONFIG.img_APRAD),
        'img_apsky': f"{CONFIG.img_APSKY[0]} {CONFIG.img_APSKY[1]}",
        'RCentroid': str(CONFIG.img_RCentroid),
        'PSFPhot': str(CONFIG.PSFPhot),
        'FSat': str(CONFIG.FSat),
        'Zero': str(CONFIG.ZERO),
    }
    
    print(f"Updating {filename} with CONFIG.py values...")
    print("=" * 60)
    
    # Read the file
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Update lines
    updated_lines = []
    changes_made = []
    
    for line in lines:
        updated = False
        original_line = line
        
        # Check if this line contains a parameter we want to update
        for param_key, new_value in updates.items():
            # Match lines like "param_key = value" or "param_key = value #comment"
            if line.strip().startswith(param_key) and '=' in line:
                # Extract comment if present
                if '#' in line:
                    comment_part = '#' + line.split('#', 1)[1]
                else:
                    comment_part = ''
                
                # Get old value
                old_value = line.split('=')[1].split('#')[0].strip()
                
                # Create new line
                new_line = f"{param_key} = {new_value}"
                
                # Pad to align comments
                if comment_part:
                    new_line = f"{new_line:40s}{comment_part}"
                else:
                    new_line = new_line + '\n'
                
                if old_value != new_value:
                    changes_made.append(f"  {param_key:15s}: {old_value:15s} → {new_value}")
                    print(f"✓ {param_key:15s}: {old_value:15s} → {new_value}")
                
                updated_lines.append(new_line)
                updated = True
                break
        
        if not updated:
            updated_lines.append(original_line)
    
    print("=" * 60)
    
    if changes_made:
        # Write back to file
        with open(filename, 'w') as f:
            f.writelines(updated_lines)
        print(f"✅ Updated {len(changes_made)} parameter(s) in {filename}")
        return True
    else:
        print(f"✓ All parameters already match CONFIG.py - no changes needed")
        return True

if __name__ == '__main__':
    success = update_dolphot_param(CONFIG.DOLPHOT_PARAM_FILE)
    sys.exit(0 if success else 1)
