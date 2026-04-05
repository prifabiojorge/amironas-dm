import re
import pprint

filename = 'e:/Documents/DEGOO/amirona/dinamica-molecular/mapeamento_projeto_20260220_152320.txt'

data = {
    'mdp': {},
    'xvg': {},
    'log': {}
}

current_section = None
current_file = None

with open(filename, 'r', encoding='utf-8', errors='replace') as f:
    for line in f:
        line = line.strip()
        
        # Detect sections
        if line.startswith('SEÇÃO 6:'):
            current_section = 6
            continue
        elif line.startswith('SEÇÃO 7:'):
            current_section = 7
            continue
        elif line.startswith('SEÇÃO 8:'):
            current_section = 8
            continue
        elif line.startswith('SEÇÃO 9:'):
            current_section = 9
            continue
            
        # Match file headers in sections
        # Example format: "--- ./alpha-amirona-6ncf/md.mdp ---" or "Arquivo: ./alpha-amirona-6ncf/md.mdp"
        file_match = re.match(r'^-{3,}\s+(.+?)\s+-{3,}$', line)
        if not file_match:
            file_match = re.match(r'^==>\s+(.+?)\s+<==$', line)
        if not file_match:
            file_match = re.search(r'(?:ARQUIVO:\s*|==>\s*|^)([\./a-zA-Z0-9_-]+\.(?:mdp|xvg|log))', line, re.IGNORECASE)
            
        if file_match and '---' in line:
             current_file = file_match.group(1).strip()
             if current_section == 6 and current_file not in data['mdp']:
                 data['mdp'][current_file] = {}
             elif current_section == 7 and current_file not in data['xvg']:
                 data['xvg'][current_file] = {'xaxis': None, 'yaxis': None}
             elif current_section == 8 and current_file not in data['log']:
                 data['log'][current_file] = {'finished': False, 'performance': None}
             continue
             
        # If no strict match but we see a filename, let's try to capture it
        if line.startswith('--- .'):
            current_file = line.replace('--- ', '').replace(' ---', '').strip()
            if current_section == 6 and current_file not in data['mdp']:
                 data['mdp'][current_file] = {}
            elif current_section == 7 and current_file not in data['xvg']:
                 data['xvg'][current_file] = {'xaxis': None, 'yaxis': None}
            elif current_section == 8 and current_file not in data['log']:
                 data['log'][current_file] = {'finished': False, 'performance': None}
            continue

        if current_section == 6 and current_file:
            # Parse mdp
            if '=' in line and not line.startswith(';'):
                parts = line.split('=')
                key = parts[0].strip()
                val_comment = '='.join(parts[1:])
                val = val_comment.split(';')[0].strip()
                if key in ['dt', 'nsteps', 'nstxout', 'nstvout', 'nstfout', 'nstxtcout', 'nstenergy', 'nstlog']:
                    data['mdp'][current_file][key] = val
                    
        elif current_section == 7 and current_file:
            # Parse xvg
            if line.startswith('@') and 'xaxis  label' in line:
                data['xvg'][current_file]['xaxis'] = line.split('"')[1] if '"' in line else line
            elif line.startswith('@') and 'yaxis  label' in line:
                data['xvg'][current_file]['yaxis'] = line.split('"')[1] if '"' in line else line
            elif line.startswith('@    xaxis  label'):
                data['xvg'][current_file]['xaxis'] = line.split('"')[1] if '"' in line else line
            elif line.startswith('@    yaxis  label'):
                data['xvg'][current_file]['yaxis'] = line.split('"')[1] if '"' in line else line
            # some xvg have exactly '@ xaxis label "Time (ps)"'
            elif re.search(r'@\s+xaxis\s+label\s+"([^"]+)"', line):
                 data['xvg'][current_file]['xaxis'] = re.search(r'@\s+xaxis\s+label\s+"([^"]+)"', line).group(1)
            elif re.search(r'@\s+yaxis\s+label\s+"([^"]+)"', line):
                 data['xvg'][current_file]['yaxis'] = re.search(r'@\s+yaxis\s+label\s+"([^"]+)"', line).group(1)

        elif current_section == 8 and current_file:
            # Parse log
            if 'Finished mdrun on' in line:
                data['log'][current_file]['finished'] = True
            elif 'Performance:' in line:
                # the next line usually has the ns/day
                data['log'][current_file]['performance'] = line

# Summarize MDP
print("=== MDP SUMMARY ===")
for f, props in data['mdp'].items():
    if 'dt' in props and 'nsteps' in props:
        try:
            dt = float(props['dt'])
            nsteps = int(props['nsteps'])
            ns = (dt * nsteps) / 1000.0  # dt is usually ps
            print(f"{f}: {ns} ns (dt={dt}, nsteps={nsteps})")
        except:
             print(f"{f}: ERROR parsing dt/nsteps ({props.get('dt')}, {props.get('nsteps')})")
    elif props:
        print(f"{f}: Missing dt or nsteps: {props}")

# Summarize XVG
print("\n=== XVG SUMMARY ===")
# limit to max 10 to avoid huge output, but group them by pattern
xvg_patterns = {}
for f, props in data['xvg'].items():
    if props['xaxis'] or props['yaxis']:
         pat = f"X: {props['xaxis']} | Y: {props['yaxis']}"
         if pat not in xvg_patterns:
             xvg_patterns[pat] = []
         xvg_patterns[pat].append(f)

for pat, files in xvg_patterns.items():
    print(f"Pattern {pat} found in {len(files)} files. Examples: {files[:2]}")

# Summarize LOG
print("\n=== LOG SUMMARY ===")
for f, props in data['log'].items():
    if props['finished']:
        print(f"{f}: FINISHED SUCCESSFULLY")
    else:
        print(f"{f}: UNKNOWN/UNFINISHED")
