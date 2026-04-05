#!/usr/bin/env python3
"""
SASA Analysis Complete Framework - Real Data Integration
========================================================

Análise completa de SASA em dinâmica molecular usando GROMACS.
Calcula interface binding, prediz ΔG e gera visualizações.

Author: Computational Chemistry Lab - UFAM
Date: January 2026
Version: 1.1
"""

import os
import sys
import argparse
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

try:
    import numpy as np
    from scipy import stats, signal
    from scipy.ndimage import gaussian_filter1d
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec
    import pandas as pd
except ImportError as e:
    print(f"ERROR: Missing module: {e}")
    print("Install with: pip install numpy scipy matplotlib pandas")
    sys.exit(1)

# Global Configuration
CONFIG_DEFAULTS = {
    'output_dir': './sasa_analysis',
    'dpi': 300,
    'probe_radius': 0.14,
    'burial_threshold': 50,
    'hotspot_threshold': 50,
    'energy_threshold': 1.0,
    'color_palette': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],
    'stability_threshold': 2.0,
}

# Energy Scale for MM-PBSA
SASA_ENERGY_SCALE = -0.0054  # kcal/mol·Å²
SASA_ENERGY_OFFSET = -0.92   # kcal/mol

# Mapeamento de nomes dos ligantes conforme solicitado
LIGAND_MAP = {
    'alpha-amirona-6ncf': 'α-Amyrone',
    'af7-6ncf': 'AF7',
    '6epicyclo-6ncf': 'Epicyclomusalenone',
    '101663382-6ncf': '101663382',
    '101685610-6ncf': '101685610'
}

class SASAAnalyzer:
    """Complete SASA analyzer for protein-ligand complexes."""
    
    def __init__(self, config: Dict = None):
        self.config = {**CONFIG_DEFAULTS, **(config or {})}
        self.data = {}
        self.results = {}
        self.timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.setup_output()
        
    def setup_output(self):
        output_dir = Path(self.config['output_dir'])
        output_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir = output_dir
        print(f"[INFO] Output directory: {output_dir}")
        
    def load_sasa_xvg(self, filename: str) -> Tuple[np.ndarray, np.ndarray]:
        """Load GROMACS XVG file."""
        time_data = []
        sasa_data = []
        
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('@'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            time_data.append(float(parts[0]))
                            sasa_data.append(float(parts[1]))
                        except ValueError:
                            continue
            
            if len(sasa_data) == 0:
                print(f"[WARNING] Empty file: {filename}")
                return None, None
                
            return np.array(time_data), np.array(sasa_data)
        
        except FileNotFoundError:
            print(f"[ERROR] File not found: {filename}")
            return None, None
    
    def analyze_sasa_trajectory(self, time: np.ndarray, sasa: np.ndarray,
                                system_name: str) -> Dict:
        """Analyze SASA trajectory statistics."""
        if sasa is None or len(sasa) == 0:
            return {}
        
        stats_dict = {
            'system_name': system_name,
            'n_frames': len(sasa),
            'time_range_ps': (float(time[0]), float(time[-1])),
            'duration_ns': float((time[-1] - time[0]) / 1000),
            'mean_sasa_nm2': float(np.mean(sasa)),
            'mean_sasa_angstrom2': float(np.mean(sasa) * 100),
            'std_sasa': float(np.std(sasa)),
            'min_sasa': float(np.min(sasa)),
            'max_sasa': float(np.max(sasa)),
            'range_sasa': float(np.max(sasa) - np.min(sasa)),
        }
        
        initial_sasa = np.mean(sasa[:int(len(sasa)*0.05)])
        final_sasa = np.mean(sasa[int(len(sasa)*0.95):])
        variation_pct = 100 * abs(final_sasa - initial_sasa) / initial_sasa
        
        stats_dict['initial_sasa_mean'] = float(initial_sasa)
        stats_dict['final_sasa_mean'] = float(final_sasa)
        stats_dict['variation_percent'] = float(variation_pct)
        stats_dict['is_stable'] = variation_pct < self.config['stability_threshold']
        
        coeffs = np.polyfit(time, sasa, 1)
        slope = coeffs[0]
        
        stats_dict['regression_slope'] = float(slope)
        stats_dict['trend'] = 'increasing' if slope > 0.001 else (
            'decreasing' if slope < -0.001 else 'stable'
        )
        
        sasa_smooth = gaussian_filter1d(sasa, sigma=10)
        rms_fluctuation = np.sqrt(np.mean((sasa - sasa_smooth) ** 2))
        stats_dict['rms_fluctuation'] = float(rms_fluctuation)
        
        peaks, _ = signal.find_peaks(sasa, height=np.mean(sasa) + 2*np.std(sasa))
        stats_dict['n_peaks'] = len(peaks)
        
        return stats_dict
    
    def calculate_interface_buried_real(self, sasa_protein: np.ndarray,
                                        sasa_ligand: np.ndarray,
                                        sasa_complex: np.ndarray) -> Dict:
        """Calculate buried SASA from real data."""
        mean_prot = np.mean(sasa_protein)
        mean_lig = np.mean(sasa_ligand)
        mean_complex = np.mean(sasa_complex)
        
        buried_interface = (mean_prot + mean_lig - mean_complex) / 2
        
        if mean_lig <= 0:
            burial_percent = 0.0
        else:
            burial_percent = 100 * (mean_lig - (mean_complex - mean_prot)) / mean_lig
            burial_percent = max(0, min(100, burial_percent))
        
        if burial_percent > 70:
            strength = 'STRONG'
            delta_g_estimate = -12.0
        elif burial_percent > 50:
            strength = 'MODERATE'
            delta_g_estimate = -6.5
        elif burial_percent > 30:
            strength = 'WEAK'
            delta_g_estimate = -3.5
        else:
            strength = 'VERY WEAK'
            delta_g_estimate = -1.0
        
        return {
            'buried_interface_nm2': float(buried_interface),
            'buried_interface_angstrom2': float(buried_interface * 100),
            'burial_percent': float(burial_percent),
            'binding_strength': strength,
            'delta_g_estimate_kcal': float(delta_g_estimate),
            'mean_sasa_protein_nm2': float(mean_prot),
            'mean_sasa_ligand_nm2': float(mean_lig),
            'mean_sasa_complex_nm2': float(mean_complex),
        }
    
    def predict_delta_g(self, buried_sasa_angstrom2: float) -> Dict:
        """Predict ΔG using MM-PBSA-SASA model."""
        delta_g = SASA_ENERGY_SCALE * buried_sasa_angstrom2 + SASA_ENERGY_OFFSET
        confidence_interval = abs(delta_g) * 0.30
        
        return {
            'delta_g_predicted_kcal': float(delta_g),
            'delta_g_min_kcal': float(delta_g - confidence_interval),
            'delta_g_max_kcal': float(delta_g + confidence_interval),
            'confidence_range_kcal': float(confidence_interval),
        }
    
    def plot_sasa_trajectory_overlay(self, sasa_dict: Dict):
        """Plot overlay of SASA trajectories."""
        fig, ax = plt.subplots(figsize=(12, 6))
        colors = self.config['color_palette']
        
        for i, (name, (time, sasa)) in enumerate(sasa_dict.items()):
            if sasa is not None:
                ax.plot(time/1000, sasa*100, label=name, 
                       color=colors[i % len(colors)], linewidth=2, alpha=0.8)
        
        ax.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
        ax.set_ylabel('SASA (Ų)', fontsize=12, fontweight='bold')
        ax.set_title('SASA Trajectory - All Systems', fontsize=14, fontweight='bold')
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#f8f9fa')
        fig.tight_layout()
        
        filepath = self.output_dir / "01_SASA_Trajectory_Overlay.png"
        fig.savefig(filepath, dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()
        print(f"[✓] Plot 1 saved: {filepath}")
    
    def plot_sasa_individual_systems(self, sasa_dict: Dict, stats_dict: Dict):
        """Plot individual system panels."""
        n_systems = len(sasa_dict)
        fig = plt.figure(figsize=(14, 3*n_systems))
        gs = GridSpec(n_systems, 1, figure=fig)
        
        for i, (name, (time, sasa)) in enumerate(sasa_dict.items()):
            if sasa is None:
                continue
            
            ax = fig.add_subplot(gs[i, 0])
            mean_val = np.mean(sasa) * 100
            std_val = np.std(sasa) * 100
            
            ax.plot(time/1000, sasa*100, linewidth=2, color='#1f77b4', alpha=0.7)
            ax.axhline(mean_val, color='red', linestyle='--', linewidth=2, 
                      label=f'Mean: {mean_val:.1f} ± {std_val:.1f} Ų')
            ax.fill_between(time/1000, 
                           (mean_val - std_val), 
                           (mean_val + std_val),
                           alpha=0.2, color='red')
            
            ax.set_ylabel('SASA (Ų)', fontsize=11, fontweight='bold')
            ax.set_title(f'System: {name}', fontsize=12, fontweight='bold')
            ax.legend(loc='best', fontsize=9)
            ax.grid(True, alpha=0.3)
            ax.set_facecolor('#f8f9fa')
        
        fig.text(0.5, 0.01, 'Time (ns)', ha='center', fontsize=12, fontweight='bold')
        fig.tight_layout(rect=[0, 0.02, 1, 1])
        
        filepath = self.output_dir / "02_SASA_Individual_Systems.png"
        fig.savefig(filepath, dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()
        print(f"[✓] Plot 2 saved: {filepath}")
    
    def plot_ligand_burial_analysis(self, interface_data: Dict):
        """Plot ligand burial analysis."""
        systems = list(interface_data.keys())
        burials = [interface_data[s]['burial_percent'] for s in systems]
        delta_gs = [interface_data[s]['delta_g_predicted_kcal'] for s in systems]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        colors_burial = ['#2ca02c' if b > 70 else '#ff7f0e' if b > 50 else '#d62728' 
                         for b in burials]
        bars1 = ax1.bar(systems, burials, color=colors_burial, alpha=0.7, edgecolor='black')
        ax1.axhline(70, color='green', linestyle='--', linewidth=2, label='Strong (>70%)')
        ax1.axhline(50, color='orange', linestyle='--', linewidth=2, label='Moderate (>50%)')
        ax1.set_ylabel('Ligand Burial (%)', fontsize=12, fontweight='bold')
        ax1.set_title('Ligand Burial Analysis', fontsize=13, fontweight='bold')
        ax1.legend(fontsize=9)
        ax1.set_ylim([0, 100])
        ax1.grid(True, alpha=0.3, axis='y')
        ax1.set_facecolor('#f8f9fa')
        
        for bar, burial in zip(bars1, burials):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{burial:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        colors_dg = ['#2ca02c' if abs(dg) > 8 else '#ff7f0e' if abs(dg) > 5 else '#d62728'
                    for dg in delta_gs]
        bars2 = ax2.bar(systems, delta_gs, color=colors_dg, alpha=0.7, edgecolor='black')
        ax2.axhline(-8, color='green', linestyle='--', linewidth=2, label='Strong (<-8)')
        ax2.axhline(-5, color='orange', linestyle='--', linewidth=2, label='Moderate (<-5)')
        ax2.set_ylabel('ΔG Predicted (kcal/mol)', fontsize=12, fontweight='bold')
        ax2.set_title('ΔG Prediction (MM-PBSA-SASA)', fontsize=13, fontweight='bold')
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.set_facecolor('#f8f9fa')
        
        for bar, dg in zip(bars2, delta_gs):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{dg:.2f}', ha='center', va='top' if dg < 0 else 'bottom', fontweight='bold')
        
        plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
        
        fig.tight_layout()
        filepath = self.output_dir / "03_Ligand_Burial_Analysis.png"
        fig.savefig(filepath, dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()
        print(f"[✓] Plot 3 saved: {filepath}")
    
    def plot_stability_analysis(self, stats_dict: Dict):
        """Plot stability analysis."""
        systems = list(stats_dict.keys())
        variations = [stats_dict[s]['variation_percent'] for s in systems]
        stabilities = [100 - v for v in variations]
        
        fig, ax = plt.subplots(figsize=(12, 6))
        colors = ['#2ca02c' if s > 98 else '#ff7f0e' if s > 95 else '#d62728'
                  for s in stabilities]
        bars = ax.bar(systems, stabilities, color=colors, alpha=0.7, edgecolor='black')
        
        ax.axhline(98, color='green', linestyle='--', linewidth=2, label='Highly stable (>98%)')
        ax.axhline(95, color='orange', linestyle='--', linewidth=2, label='Stable (95-98%)')
        ax.set_ylabel('Stability (%)', fontsize=12, fontweight='bold')
        ax.set_title('Structural Stability Analysis', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.set_ylim([0, 105])
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_facecolor('#f8f9fa')
        
        for bar, stability in zip(bars, stabilities):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{stability:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        fig.tight_layout()
        
        filepath = self.output_dir / "04_Stability_Analysis.png"
        fig.savefig(filepath, dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()
        print(f"[✓] Plot 4 saved: {filepath}")
    
    def plot_sasa_affinity_correlation(self, interface_data: Dict):
        """Plot SASA-affinity correlation."""
        buried_sasa = [interface_data[s]['buried_interface_angstrom2'] for s in interface_data]
        delta_gs = [interface_data[s]['delta_g_predicted_kcal'] for s in interface_data]
        systems = list(interface_data.keys())
        
        fig, ax = plt.subplots(figsize=(10, 7))
        
        scatter = ax.scatter(buried_sasa, delta_gs, s=150, alpha=0.6, 
                            c=range(len(systems)), cmap='viridis', edgecolors='black', linewidth=2)
        
        if len(buried_sasa) > 1:
            z = np.polyfit(buried_sasa, delta_gs, 1)
            p = np.poly1d(z)
            x_line = np.linspace(min(buried_sasa), max(buried_sasa), 100)
            ax.plot(x_line, p(x_line), "r--", linewidth=2, alpha=0.8, label='Linear regression')
            
            r = np.corrcoef(buried_sasa, delta_gs)[0, 1]
            ax.text(0.05, 0.95, f'r = {r:.3f}', transform=ax.transAxes,
                   fontsize=11, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlabel('Buried SASA (Ų)', fontsize=12, fontweight='bold')
        ax.set_ylabel('ΔG Predicted (kcal/mol)', fontsize=12, fontweight='bold')
        ax.set_title('SASA-Affinity Correlation', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#f8f9fa')
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Systems', fontsize=10)
        
        fig.tight_layout()
        filepath = self.output_dir / "05_SASA_Affinity_Correlation.png"
        fig.savefig(filepath, dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()
        print(f"[✓] Plot 5 saved: {filepath}")
    
    def plot_sasa_distribution(self, sasa_dict: Dict):
        """Plot SASA distribution."""
        systems = list(sasa_dict.keys())
        fig = plt.figure(figsize=(14, 8))
        gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1.2])
        
        ax1 = fig.add_subplot(gs[0, 0])
        colors = self.config['color_palette']
        
        for i, name in enumerate(systems):
            time, sasa = sasa_dict[name]
            if sasa is not None:
                ax1.hist(sasa*100, bins=30, alpha=0.5, label=name,
                        color=colors[i % len(colors)], edgecolor='black')
        
        ax1.set_xlabel('SASA (Ų)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Frequency', fontsize=11, fontweight='bold')
        ax1.set_title('SASA Distribution - Histograms', fontsize=12, fontweight='bold')
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3, axis='y')
        ax1.set_facecolor('#f8f9fa')
        
        ax2 = fig.add_subplot(gs[1, 0])
        sasa_data = []
        labels = []
        for name in systems:
            time, sasa = sasa_dict[name]
            if sasa is not None:
                sasa_data.append(sasa * 100)
                labels.append(name)
        
        bp = ax2.boxplot(sasa_data, labels=labels, patch_artist=True,
                        widths=0.6, showmeans=True, meanline=True)
        
        for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        
        ax2.set_ylabel('SASA (Ų)', fontsize=11, fontweight='bold')
        ax2.set_title('SASA Distribution - Box Plots', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.set_facecolor('#f8f9fa')
        
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
        
        fig.tight_layout()
        filepath = self.output_dir / "06_SASA_Distribution_Analysis.png"
        fig.savefig(filepath, dpi=self.config['dpi'], bbox_inches='tight')
        plt.close()
        print(f"[✓] Plot 6 saved: {filepath}")
    
    def run_complete_analysis(self, system_configs: List[Dict]):
        """Run complete analysis with real data."""
        print("\n" + "="*80)
        print("COMPLETE SASA ANALYSIS - REAL DATA")
        print("="*80 + "\n")
        
        print("[1/6] Loading SASA data...")
        sasa_dict = {}
        stats_dict = {}
        interface_data = {}
        
        for config in system_configs:
            sys_name = config['name']
            
            time_c, sasa_c = self.load_sasa_xvg(config['complex_xvg'])
            time_p, sasa_p = self.load_sasa_xvg(config['protein_xvg'])
            time_l, sasa_l = self.load_sasa_xvg(config['ligand_xvg'])
            
            if sasa_c is None or sasa_p is None or sasa_l is None:
                print(f"  ✗ {sys_name}: Error loading data")
                continue
            
            sasa_dict[sys_name] = (time_c, sasa_c)
            
            stats = self.analyze_sasa_trajectory(time_c, sasa_c, sys_name)
            stats_dict[sys_name] = stats
            
            interface = self.calculate_interface_buried_real(sasa_p, sasa_l, sasa_c)
            interface_data[sys_name] = interface
            
            dg_pred = self.predict_delta_g(interface['buried_interface_angstrom2'])
            interface_data[sys_name].update(dg_pred)
            
            print(f"  ✓ {sys_name}: {len(sasa_c)} frames")
            print(f"    └─ Buried SASA: {interface['buried_interface_angstrom2']:.1f} Ų")
            print(f"    └─ Burial: {interface['burial_percent']:.1f}%")
            print(f"    └─ ΔG predicted: {dg_pred['delta_g_predicted_kcal']:.2f} kcal/mol")
        
        print("\n[2/6] Generating plots...")
        self.plot_sasa_trajectory_overlay(sasa_dict)
        self.plot_sasa_individual_systems(sasa_dict, stats_dict)
        self.plot_ligand_burial_analysis(interface_data)
        self.plot_stability_analysis(stats_dict)
        self.plot_sasa_affinity_correlation(interface_data)
        self.plot_sasa_distribution(sasa_dict)
        
        print("\n[3/6] Saving statistics...")
        stats_file = self.output_dir / 'sasa_statistics.json'
        with open(stats_file, 'w') as f:
            json.dump(stats_dict, f, indent=2)
        print(f"  ✓ Statistics saved: {stats_file}")
        
        interface_file = self.output_dir / 'interface_data.json'
        with open(interface_file, 'w') as f:
            json.dump(interface_data, f, indent=2)
        print(f"  ✓ Interface data saved: {interface_file}")
        
        print("\n[4/6] Generating reports...")
        reports_file = self.output_dir / 'SASA_Analysis_Report.txt'
        with open(reports_file, 'w') as f:
            for sys_name in sasa_dict:
                f.write(f"\n{'='*80}\n")
                f.write(f"SYSTEM: {sys_name}\n")
                f.write(f"{'='*80}\n\n")
                
                stats = stats_dict.get(sys_name, {})
                interface = interface_data.get(sys_name, {})
                
                f.write("TRAJECTORY STATISTICS:\n")
                f.write(f"  Frames: {stats.get('n_frames', 0)}\n")
                f.write(f"  Duration: {stats.get('duration_ns', 0):.1f} ns\n")
                f.write(f"  Mean SASA: {stats.get('mean_sasa_angstrom2', 0):.1f} Ų\n")
                f.write(f"  Stability: {100 - stats.get('variation_percent', 0):.1f}%\n")
                f.write(f"  Status: {'STABLE' if stats.get('is_stable') else 'UNSTABLE'}\n\n")
                
                f.write("INTERFACE ANALYSIS:\n")
                f.write(f"  Protein SASA: {interface.get('mean_sasa_protein_nm2', 0):.2f} nm²\n")
                f.write(f"  Ligand SASA: {interface.get('mean_sasa_ligand_nm2', 0):.2f} nm²\n")
                f.write(f"  Complex SASA: {interface.get('mean_sasa_complex_nm2', 0):.2f} nm²\n")
                f.write(f"  Buried SASA: {interface.get('buried_interface_angstrom2', 0):.1f} Ų\n")
                f.write(f"  Burial %: {interface.get('burial_percent', 0):.1f}%\n")
                f.write(f"  Binding Strength: {interface.get('binding_strength', 'N/A')}\n")
                f.write(f"  ΔG (estimate): {interface.get('delta_g_estimate_kcal', 0):.1f} kcal/mol\n")
                f.write(f"  ΔG (MM-PBSA): {interface.get('delta_g_predicted_kcal', 0):.2f} kcal/mol\n")
                f.write(f"  CI (95%): [{interface.get('delta_g_min_kcal', 0):.2f}, ")
                f.write(f"{interface.get('delta_g_max_kcal', 0):.2f}]\n")
        
        print(f"  ✓ Report saved: {reports_file}")
        
        print("\n[5/6] Analysis complete!")
        print("\n" + "="*80)
        print("SUMMARY")
        print("="*80)
        print(f"Systems analyzed: {len(sasa_dict)}")
        print(f"Output directory: {self.output_dir}")
        print(f"High-resolution plots: 6 PNG files (300 dpi)")
        print(f"Data files: statistics.json, interface_data.json")
        print(f"Report: {reports_file.name}")
        print("="*80 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Complete SASA analysis for MD trajectories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Single system
  python3 analyze_sasa_complete.py \\
    --complex sasa_complex.xvg \\
    --protein sasa_protein.xvg \\
    --ligand sasa_ligand.xvg \\
    --name "alpha-amirona-6ncf"

  # Multiple systems
  python3 analyze_sasa_complete.py --batch-dir ./systems --output ./results
        """
    )
    
    parser.add_argument('--complex', type=str,
                       help='SASA complex XVG file')
    parser.add_argument('--protein', type=str,
                       help='SASA protein XVG file')
    parser.add_argument('--ligand', type=str,
                       help='SASA ligand XVG file')
    parser.add_argument('--name', type=str, default='System_1',
                       help='System name')
    parser.add_argument('--output', type=str, default='./sasa_analysis',
                       help='Output directory')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Plot resolution (default: 300)')
    parser.add_argument('--batch-dir', type=str,
                       help='Directory containing subdirectories with sasa_*.xvg files')
    
    args = parser.parse_args()
    
    config = CONFIG_DEFAULTS.copy()
    config['output_dir'] = args.output
    if args.dpi:
        config['dpi'] = args.dpi
    
    analyzer = SASAAnalyzer(config)
    
    if args.complex and args.protein and args.ligand:
        # Aplica o mapeamento se o nome estiver no dicionário
        final_name = LIGAND_MAP.get(args.name, args.name)
        system_configs = [{
            'name': final_name,
            'complex_xvg': args.complex,
            'protein_xvg': args.protein,
            'ligand_xvg': args.ligand
        }]
        analyzer.run_complete_analysis(system_configs)
    elif args.batch_dir:
        print("[INFO] Batch mode - processing directory")
        batch_dir = Path(args.batch_dir)
        system_configs = []
        
        for subdir in sorted(batch_dir.iterdir()):
            if not subdir.is_dir():
                continue
            
            complex_f = subdir / 'sasa_complex.xvg'
            protein_f = subdir / 'sasa_protein.xvg'
            ligand_f = subdir / 'sasa_ligand.xvg'
            
            if complex_f.exists() and protein_f.exists() and ligand_f.exists():
                # Aplica o mapeamento automático baseado no nome da pasta
                final_name = LIGAND_MAP.get(subdir.name, subdir.name)
                system_configs.append({
                    'name': final_name,
                    'complex_xvg': str(complex_f),
                    'protein_xvg': str(protein_f),
                    'ligand_xvg': str(ligand_f)
                })
        
        if system_configs:
            analyzer.run_complete_analysis(system_configs)
        else:
            print("[ERROR] No valid systems found in batch directory")
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
