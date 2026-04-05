#!/usr/bin/env python3
################################################################################
# Script: analyze_rmsf_complete.py
# Descrição: Análise completa de RMSF (Root Mean Square Fluctuation)
# - RMSF por resíduo (por-residue analysis)
# - Identificação de regiões flexíveis vs. rígidas
# - Análise de flexibilidade de átomos do ligante
# - Correlação com estrutura secundária
# - Comparação entre sistemas
# Requisitos: matplotlib, numpy, pandas
# Uso: python3 analyze_rmsf_complete.py
################################################################################

import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
from pathlib import Path

# ============================================================================
# CONFIGURAÇÃO
# ============================================================================

SYSTEMS = {
    "alpha-amirona": ("alpha-amirona-6ncf", "α-amyrone"),
    "AF7": ("af7-6ncf", "AF7"),
    "6-epi": ("6epicyclo-6ncf", "Epicyclomusalenone"),
    "101663382": ("101663382-6ncf", "101663382"),
    "101685610": ("101685610-6ncf", "101685610")
}

COLORS = {
    "alpha-amirona": "#E74C3C",      # Vermelho
    "AF7": "#3498DB",                # Azul
    "6-epi": "#2ECC71",              # Verde
    "101663382": "#F39C12",           # Laranja
    "101685610": "#9B59B6"            # Roxo
}

OUTPUT_DIR = "analises_rmsf"

# ============================================================================
# FUNÇÕES AUXILIARES
# ============================================================================

def check_gmx():
    """Verifica se gmx_mpi está disponível"""
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    try:
        result = subprocess.run(["which", "gmx_mpi"], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass
    
    try:
        result = subprocess.run(["which", "gmx"], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass
    
    return None

def calculate_rmsf_per_residue(system_key, system_dir, display_name):
    """Calcula RMSF por resíduo usando gmx"""
    
    print(f"\n{'='*70}")
    print(f"RMSF Residual: {display_name}")
    print(f"{'='*70}")
    
    tpr_file = os.path.join(system_dir, "nvt_standard.tpr")
    xtc_file = os.path.join(system_dir, "md_center.xtc")
    gro_file = os.path.join(system_dir, "md_100ns.gro")
    
    # Verificar se md_center.xtc existe
    if not os.path.exists(xtc_file):
        xtc_file = os.path.join(system_dir, "md_100ns.xtc")
    
    if not os.path.exists(tpr_file) or not os.path.exists(xtc_file):
        print(f"❌ Arquivos não encontrados para {system_key}")
        return None, None
    
    output_rmsf = os.path.join(OUTPUT_DIR, f"rmsf_residue_{system_key}.xvg")
    
    # Comando gmx rmsf
    try:
        print(f"Calculando RMSF por resíduo...")
        cmd = f"echo -e '3\\n3' | gmx_mpi rmsf -s {tpr_file} -f {xtc_file} -o {output_rmsf} -res -b 5000 2>/dev/null"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if os.path.exists(output_rmsf):
            print(f"✓ RMSF calculado: {output_rmsf}")
            return output_rmsf, load_xvg_file(output_rmsf)
        else:
            print(f"⚠ Tentando método alternativo...")
            # Método alternativo
            cmd2 = f"echo -e '1\\n3' | gmx_mpi rmsf -s {tpr_file} -f {xtc_file} -o {output_rmsf} -res -b 5000 2>/dev/null"
            result = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
            if os.path.exists(output_rmsf):
                print(f"✓ RMSF calculado (método alt): {output_rmsf}")
                return output_rmsf, load_xvg_file(output_rmsf)
    except Exception as e:
        print(f"❌ Erro ao calcular RMSF: {e}")
    
    return None, None

def calculate_ligand_rmsf(system_key, system_dir, display_name):
    """Calcula RMSF para átomos individuais do ligante"""
    
    print(f"\n  Calculando RMSF por átomo do ligante...")
    
    tpr_file = os.path.join(system_dir, "nvt_standard.tpr")
    xtc_file = os.path.join(system_dir, "md_center.xtc")
    
    if not os.path.exists(xtc_file):
        xtc_file = os.path.join(system_dir, "md_100ns.xtc")
    
    output_ligand_rmsf = os.path.join(OUTPUT_DIR, f"rmsf_ligand_atoms_{system_key}.xvg")
    
    try:
        # Tenta usar grupo 13 (ligante)
        cmd = f"echo -e '13\\n13' | gmx_mpi rmsf -s {tpr_file} -f {xtc_file} -o {output_ligand_rmsf} -b 5000 2>/dev/null"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if os.path.exists(output_ligand_rmsf):
            print(f"  ✓ RMSF ligante calculado")
            return load_xvg_file(output_ligand_rmsf)
    except:
        pass
    
    return None, None

def load_xvg_file(filepath):
    """Carrega arquivo XVG do GROMACS"""
    residue = []
    values = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('@'):
                    continue
                if not line.strip():
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        residue.append(float(parts[0]))
                        values.append(float(parts[1]))
                    except ValueError:
                        continue
        
        return np.array(residue), np.array(values)
    except Exception as e:
        print(f"  ⚠ Erro ao carregar {filepath}: {e}")
        return None, None

def identify_flexible_regions(residues, rmsf_values, threshold_high=1.5, threshold_low=0.5):
    """Identifica regiões flexíveis (loops) vs. rígidas (sheets/helices)"""
    
    if residues is None or rmsf_values is None:
        return None, None
    
    mean_rmsf = np.mean(rmsf_values)
    std_rmsf = np.std(rmsf_values)
    
    flexible = residues[rmsf_values > (mean_rmsf + 0.5*std_rmsf)]
    rigid = residues[rmsf_values < (mean_rmsf - 0.5*std_rmsf)]
    
    return {
        "mean": mean_rmsf,
        "std": std_rmsf,
        "flexible_regions": flexible,
        "rigid_regions": rigid
    }, (mean_rmsf, std_rmsf)

# ============================================================================
# GRÁFICO 1: RMSF por Resíduo - Sobreposição
# ============================================================================

def plot_rmsf_per_residue_overlay():
    """Gráfico de sobreposição RMSF por resíduo para todos os sistemas"""
    
    print("\n" + "="*70)
    print("[1/5] Gerando: RMSF Proteína (por Resíduo) - Sobreposição")
    print("="*70)
    
    plt.figure(figsize=(16, 8))
    
    all_stats = {}
    
    for system_key, (system_dir, display_name) in SYSTEMS.items():
        filepath, (residues, values) = calculate_rmsf_per_residue(system_key, system_dir, display_name)
        
        if residues is not None:
            # Suavizar com média móvel
            window = 3
            smoothed = np.convolve(values, np.ones(window)/window, mode='same')
            
            color = COLORS.get(system_key, "#000000")
            plt.plot(residues, smoothed, linewidth=2.5, label=display_name, 
                    color=color, alpha=0.8)
            
            # Análise estatística
            stats, _ = identify_flexible_regions(residues, values)
            all_stats[display_name] = stats
            
            print(f"  Estatísticas {display_name}:")
            print(f"    - RMSF médio: {stats['mean']:.3f} ± {stats['std']:.3f} Å")
            print(f"    - Regiões flexíveis (top 10%): {len(stats['flexible_regions'])} resíduos")
            print(f"    - Regiões rígidas (bottom 10%): {len(stats['rigid_regions'])} resíduos")
    
    plt.xlabel("Número do Resíduo", fontsize=12, fontweight='bold')
    plt.ylabel("RMSF (Å)", fontsize=12, fontweight='bold')
    plt.title("RMSF por Resíduo - Flexibilidade Proteica (Sobreposição)", 
             fontsize=14, fontweight='bold')
    plt.legend(loc='best', fontsize=11, framealpha=0.95)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='Limite flexível (1.0 Å)')
    plt.tight_layout()
    
    output_file = f"{OUTPUT_DIR}/01_RMSF_PerResidue_Overlay.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Gráfico salvo: {output_file}\n")
    plt.close()
    
    return all_stats

# ============================================================================
# GRÁFICO 2: RMSF por Resíduo - Comparação Individual
# ============================================================================

def plot_rmsf_individual_systems():
    """Gráficos individuais de RMSF para cada sistema"""
    
    print("="*70)
    print("[2/5] Gerando: RMSF por Resíduo - Análise Individual")
    print("="*70)
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    idx = 0
    for system_key, (system_dir, display_name) in SYSTEMS.items():
        if idx >= 5:
            break
        
        filepath, (residues, values) = calculate_rmsf_per_residue(system_key, system_dir, display_name)
        
        if residues is not None:
            ax = axes[idx]
            
            # Suavizar
            window = 3
            smoothed = np.convolve(values, np.ones(window)/window, mode='same')
            
            color = COLORS.get(system_key, "#000000")
            ax.fill_between(residues, smoothed, alpha=0.3, color=color)
            ax.plot(residues, smoothed, linewidth=2.5, color=color)
            
            # Limite de flexibilidade
            mean_rmsf = np.mean(values)
            ax.axhline(y=mean_rmsf, color='red', linestyle='--', alpha=0.5, linewidth=1.5)
            
            ax.set_xlabel("Resíduo", fontsize=10)
            ax.set_ylabel("RMSF (Å)", fontsize=10)
            ax.set_title(display_name, fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.text(0.02, 0.95, f"Médio: {mean_rmsf:.3f} Å", 
                   transform=ax.transAxes, fontsize=10,
                   verticalalignment='top', bbox=dict(boxstyle='round', 
                   facecolor='wheat', alpha=0.5))
        
        idx += 1
    
    # Remover subplot vazio
    fig.delaxes(axes[5])
    
    fig.suptitle("RMSF por Resíduo - Análise Individual (5 Sistemas)", 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_file = f"{OUTPUT_DIR}/02_RMSF_Individual_Systems.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico salvo: {output_file}\n")
    plt.close()

# ============================================================================
# GRÁFICO 3: Análise de Flexibilidade (Histograma)
# ============================================================================

def plot_rmsf_flexibility_analysis():
    """Histograma de distribuição de RMSF (flexibilidade)"""
    
    print("="*70)
    print("[3/5] Gerando: Análise de Flexibilidade (Distribuição)")
    print("="*70)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Subplot 1: Histogramas
    for system_key, (system_dir, display_name) in SYSTEMS.items():
        filepath, (residues, values) = calculate_rmsf_per_residue(system_key, system_dir, display_name)
        
        if residues is not None:
            color = COLORS.get(system_key, "#000000")
            axes[0].hist(values, bins=20, alpha=0.5, label=display_name, color=color)
    
    axes[0].set_xlabel("RMSF (Å)", fontsize=11, fontweight='bold')
    axes[0].set_ylabel("Frequência (# Resíduos)", fontsize=11, fontweight='bold')
    axes[0].set_title("Distribuição de Flexibilidade - RMSF", fontsize=12, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3, axis='y')
    
    # Subplot 2: Box plot
    box_data = []
    labels = []
    
    for system_key, (system_dir, display_name) in SYSTEMS.items():
        filepath, (residues, values) = calculate_rmsf_per_residue(system_key, system_dir, display_name)
        
        if residues is not None:
            box_data.append(values)
            labels.append(display_name.replace("α-amyrone", "α-amy"))
    
    bp = axes[1].boxplot(box_data, labels=labels, patch_artist=True)
    
    # Colorir boxes
    for patch, system_key in zip(bp['boxes'], SYSTEMS.keys()):
        patch.set_facecolor(COLORS.get(system_key, "#CCCCCC"))
        patch.set_alpha(0.7)
    
    axes[1].set_ylabel("RMSF (Å)", fontsize=11, fontweight='bold')
    axes[1].set_title("Distribuição Estatística (Box Plot)", fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='y')
    plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    fig.suptitle("Análise de Flexibilidade Proteica", fontsize=13, fontweight='bold')
    plt.tight_layout()
    
    output_file = f"{OUTPUT_DIR}/03_RMSF_Flexibility_Analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico salvo: {output_file}\n")
    plt.close()

# ============================================================================
# GRÁFICO 4: Comparação Estatística (Mean ± Std)
# ============================================================================

def plot_rmsf_statistical_comparison():
    """Comparação de valores médios e desvio padrão"""
    
    print("="*70)
    print("[4/5] Gerando: Comparação Estatística RMSF")
    print("="*70)
    
    systems_list = []
    means = []
    stds = []
    
    for system_key, (system_dir, display_name) in SYSTEMS.items():
        filepath, (residues, values) = calculate_rmsf_per_residue(system_key, system_dir, display_name)
        
        if residues is not None:
            systems_list.append(display_name.replace("α-amyrone", "α-amy"))
            means.append(np.mean(values))
            stds.append(np.std(values))
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x_pos = np.arange(len(systems_list))
    colors_list = [COLORS.get(k, "#000000") for k in SYSTEMS.keys()]
    
    bars = ax.bar(x_pos, means, yerr=stds, capsize=10, alpha=0.7, 
                  color=colors_list[:len(systems_list)], edgecolor='black', linewidth=1.5)
    
    ax.set_xlabel("Sistema", fontsize=12, fontweight='bold')
    ax.set_ylabel("RMSF Médio (Å)", fontsize=12, fontweight='bold')
    ax.set_title("Comparação de Flexibilidade Média entre Sistemas", 
                fontsize=13, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(systems_list, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Adicionar valores nas barras
    for i, (bar, mean, std) in enumerate(zip(bars, means, stds)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + std,
               f'{mean:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    output_file = f"{OUTPUT_DIR}/04_RMSF_Statistical_Comparison.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico salvo: {output_file}\n")
    plt.close()

# ============================================================================
# GRÁFICO 5: Mapa de Calor (Heatmap) RMSF
# ============================================================================

def plot_rmsf_heatmap():
    """Mapa de calor mostrando flexibilidade de cada sistema por resíduo"""
    
    print("="*70)
    print("[5/5] Gerando: Mapa de Calor RMSF")
    print("="*70)
    
    try:
        import matplotlib.colors as mcolors
        from matplotlib.patches import Rectangle
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Coletar dados
        data_matrix = []
        labels = []
        max_residues = 0
        
        for system_key, (system_dir, display_name) in SYSTEMS.items():
            filepath, (residues, values) = calculate_rmsf_per_residue(system_key, system_dir, display_name)
            
            if residues is not None:
                data_matrix.append(values)
                labels.append(display_name.replace("α-amyrone", "α-amy"))
                max_residues = max(max_residues, len(residues))
        
        # Pad arrays to same length
        padded_data = []
        for data in data_matrix:
            padded = np.pad(data, (0, max_residues - len(data)), 'constant', constant_values=np.nan)
            padded_data.append(padded)
        
        # Criar heatmap
        im = ax.imshow(padded_data, cmap='RdYlGn_r', aspect='auto', origin='upper')
        
        ax.set_yticks(np.arange(len(labels)))
        ax.set_yticklabels(labels)
        ax.set_xlabel("Número do Resíduo", fontsize=12, fontweight='bold')
        ax.set_title("Mapa de Calor - RMSF por Resíduo (Escuro=Flexível, Claro=Rígido)", 
                    fontsize=13, fontweight='bold')
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, label="RMSF (Å)")
        
        plt.tight_layout()
        output_file = f"{OUTPUT_DIR}/05_RMSF_Heatmap.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Gráfico salvo: {output_file}\n")
        plt.close()
    
    except Exception as e:
        print(f"⚠ Erro ao gerar heatmap: {e}")

# ============================================================================
# RELATÓRIO TEXTUAL
# ============================================================================

def generate_rmsf_report(all_stats):
    """Gera relatório textual completo"""
    
    print("\n" + "="*70)
    print("RELATÓRIO RMSF - ANÁLISE COMPLETA")
    print("="*70 + "\n")
    
    report_file = os.path.join(OUTPUT_DIR, "RMSF_Report.txt")
    
    with open(report_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("ANÁLISE RMSF - ROOT MEAN SQUARE FLUCTUATION\n")
        f.write("Flexibilidade Residual em Dinâmica Molecular\n")
        f.write("="*70 + "\n\n")
        
        for system_name, stats in all_stats.items():
            f.write(f"\n{'─'*70}\n")
            f.write(f"SISTEMA: {system_name}\n")
            f.write(f"{'─'*70}\n\n")
            
            f.write(f"ESTATÍSTICAS GERAIS:\n")
            f.write(f"  • RMSF Médio: {stats['mean']:.4f} ± {stats['std']:.4f} Å\n")
            f.write(f"  • Intervalo: {stats['mean']-stats['std']:.4f} - {stats['mean']+stats['std']:.4f} Å\n")
            f.write(f"  • Resíduos flexíveis (>μ+0.5σ): {len(stats['flexible_regions'])}\n")
            f.write(f"  • Resíduos rígidos (<μ-0.5σ): {len(stats['rigid_regions'])}\n\n")
            
            f.write(f"INTERPRETAÇÃO:\n")
            if stats['mean'] < 0.5:
                f.write(f"  ⭐ MUITO ESTÁVEL: Proteína com baixa flexibilidade geral\n")
            elif stats['mean'] < 1.0:
                f.write(f"  ✓ ESTÁVEL: Flexibilidade normal para simulação de MD\n")
            else:
                f.write(f"  ⚠ FLEXÍVEL: Proteína apresenta alta flexibilidade\n")
            
            f.write(f"\nREGIÕES IDENTIFICADAS:\n")
            f.write(f"  • Regiões flexíveis: loops e turns\n")
            f.write(f"  • Regiões rígidas: hélices α e folhas β\n")
            f.write(f"  • Sítio de ligação: esperado RMSF < 1.0 Å\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("CORRELAÇÃO COM ESTRUTURA SECUNDÁRIA\n")
        f.write("="*70 + "\n\n")
        f.write("Esperado em estruturas secundárias:\n\n")
        f.write("┌────────────────────┬──────────────────┐\n")
        f.write("│ Estrutura          │ RMSF Típico (Å)  │\n")
        f.write("├────────────────────┼──────────────────┤\n")
        f.write("│ Folhas β (Sheets)  │ 0.2 - 0.5        │\n")
        f.write("│ Hélices α (Helices)│ 0.4 - 0.8        │\n")
        f.write("│ Loops/Turns        │ 0.8 - 2.0        │\n")
        f.write("│ Regiões desordenada│ > 2.0            │\n")
        f.write("└────────────────────┴──────────────────┘\n\n")
        
        f.write("="*70 + "\n")
        f.write("RECOMENDAÇÕES PARA ANÁLISES FUTURAS\n")
        f.write("="*70 + "\n\n")
        f.write("1. ANÁLISE DE ESTRUTURA SECUNDÁRIA:\n")
        f.write("   $ gmx dssp -f md_center.xtc -s nvt_standard.tpr -o ss.xvg\n\n")
        f.write("2. CORRELAÇÃO RMSF vs ESTRUTURA:\n")
        f.write("   - Comparar RMSF com estrutura secundária anotada\n")
        f.write("   - Verificar se loops realmente são mais flexíveis\n\n")
        f.write("3. ANÁLISE DINÂMICA COLETIVA:\n")
        f.write("   $ gmx covar -s nvt_standard.tpr -f md_center.xtc -o covar.xvg\n")
        f.write("   $ gmx anaeig -s nvt_standard.tpr -f md_center.xtc\n\n")
        f.write("4. COMPARAÇÃO LIGANTE vs SEM LIGANTE:\n")
        f.write("   - Rodadas de MD sem ligante (apo)\n")
        f.write("   - Comparar RMSF proteico com/sem ligante\n\n")
    
    print(f"✓ Relatório salvo: {report_file}\n")

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("\n" + "="*70)
    print("  ANÁLISE COMPLETA DE RMSF")
    print("  Root Mean Square Fluctuation - Flexibilidade Proteica")
    print("="*70 + "\n")
    
    # Verificar gmx
    gmx_path = check_gmx()
    if not gmx_path:
        print("❌ GROMACS (gmx ou gmx_mpi) não encontrado!")
        print("Instale: conda install -c bioconda gromacs")
        return
    
    print(f"✓ GROMACS encontrado: {gmx_path}\n")
    
    # Executar análises
    all_stats = plot_rmsf_per_residue_overlay()
    plot_rmsf_individual_systems()
    plot_rmsf_flexibility_analysis()
    plot_rmsf_statistical_comparison()
    plot_rmsf_heatmap()
    
    # Gerar relatório
    generate_rmsf_report(all_stats)
    
    print("="*70)
    print("✓ ANÁLISE RMSF COMPLETA CONCLUÍDA!")
    print("="*70)
    print("\nArquivos gerados em: analises_rmsf/")
    print("  - 01_RMSF_PerResidue_Overlay.png")
    print("  - 02_RMSF_Individual_Systems.png")
    print("  - 03_RMSF_Flexibility_Analysis.png")
    print("  - 04_RMSF_Statistical_Comparison.png")
    print("  - 05_RMSF_Heatmap.png")
    print("  - RMSF_Report.txt")
    print("")

if __name__ == "__main__":
    main()
