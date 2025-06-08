from shiny import App, render, ui, reactive
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve
import math

# Acetone-Ethanol VLE data (approximate)
def equilibrium_curve(x):
    """
    Acetone-Ethanol equilibrium relationship
    Based on relative volatility approximation
    """
    alpha = 2.4  # Relative volatility of acetone to ethanol
    return (alpha * x) / (1 + (alpha - 1) * x)

def x_from_y_equilibrium(y):
    """Inverse equilibrium relationship"""
    alpha = 2.4
    return y / (alpha - y * (alpha - 1))

app_ui = ui.page_fluid(
    ui.h2("McCabe-Thiele Distillation Design Tool"),
    ui.h4("Acetone-Ethanol Separation"),
    ui.br(),
    
    ui.row(
        ui.column(4,
            ui.card(
                ui.card_header("Operating Conditions"),
                ui.input_numeric("feed_rate", "Feed Rate (kmol/hr):", value=100, min=10, max=1000),
                ui.input_numeric("x_feed", "Feed Composition (mol fraction acetone):", value=0.5, min=0.01, max=0.99, step=0.01),
                ui.input_numeric("x_distillate", "Distillate Composition:", value=0.95, min=0.5, max=0.99, step=0.01),
                ui.input_numeric("x_bottoms", "Bottoms Composition:", value=0.05, min=0.01, max=0.5, step=0.01),
                ui.input_numeric("q_factor", "q-factor (mol vapor condensed/mol feed):", value=1.167, min=0.5, max=2.0, step=0.01),
                ui.input_numeric("reflux_ratio_factor", "Reflux Ratio Factor (R/Rmin):", value=1.25, min=1.1, max=3.0, step=0.05),
            )
        ),
        
        ui.column(8,
            ui.output_plot("mccabe_thiele_plot", height="600px")
        )
    ),
    
    ui.br(),
    
    ui.row(
        ui.column(6,
            ui.card(
                ui.card_header("Design Results"),
                ui.output_text_verbatim("design_results")
            )
        ),
        ui.column(6,
            ui.card(
                ui.card_header("Material Balance"),
                ui.output_text_verbatim("material_balance")
            )
        )
    ),
    
    ui.br(),
    
    ui.row(
        ui.column(12,
            ui.card(
                ui.card_header("Stage-by-Stage Analysis"),
                ui.output_table("stage_analysis")
            )
        )
    )
)

def server(input, output, session):
    
    @reactive.calc
    def calculate_design():
        # Input parameters
        F = input.feed_rate()
        x_f = input.x_feed()
        x_d = input.x_distillate()
        x_b = input.x_bottoms()
        q = input.q_factor()
        R_factor = input.reflux_ratio_factor()
        
        # Overall material balance
        # F = D + B
        # F * x_f = D * x_d + B * x_b
        D = F * (x_f - x_b) / (x_d - x_b)
        B = F - D
        
        # Find minimum reflux ratio
        # At minimum reflux, operating lines intersect at feed line intersection with equilibrium curve
        
        # Feed line: y = q/(q-1) * x - x_f/(q-1)
        # For subcooled liquid, q > 1
        
        # Find intersection of feed line with equilibrium curve
        def feed_eq_intersection(x):
            y_feed = q/(q-1) * x - x_f/(q-1)
            y_eq = equilibrium_curve(x)
            return y_feed - y_eq
        
        try:
            x_intersect = fsolve(feed_eq_intersection, x_f)[0]
            y_intersect = equilibrium_curve(x_intersect)
        except:
            x_intersect = x_f
            y_intersect = equilibrium_curve(x_f)
        
        # Minimum reflux ratio
        R_min = (x_d - y_intersect) / (y_intersect - x_intersect)
        if R_min < 0:
            R_min = 0.5
        
        # Actual reflux ratio
        R = R_factor * R_min
        
        # Operating line slopes and intercepts
        # Rectifying section: y = (R/(R+1)) * x + x_d/(R+1)
        slope_rect = R / (R + 1)
        intercept_rect = x_d / (R + 1)
        
        # Stripping section operating line intersects rectifying line at feed stage
        # Find intersection of rectifying line with feed line
        # Rectifying: y = slope_rect * x + intercept_rect
        # Feed: y = q/(q-1) * x - x_f/(q-1)
        
        def find_feed_intersection(x):
            y_rect = slope_rect * x + intercept_rect
            y_feed = q/(q-1) * x - x_f/(q-1)
            return y_rect - y_feed
        
        try:
            x_feed_intersect = fsolve(find_feed_intersection, x_f)[0]
            y_feed_intersect = slope_rect * x_feed_intersect + intercept_rect
        except:
            x_feed_intersect = x_f
            y_feed_intersect = slope_rect * x_f + intercept_rect
        
        # Stripping section: passes through (x_b, x_b) and (x_feed_intersect, y_feed_intersect)
        if x_feed_intersect != x_b:
            slope_strip = (y_feed_intersect - x_b) / (x_feed_intersect - x_b)
        else:
            slope_strip = 1.0
        intercept_strip = x_b - slope_strip * x_b
        
        # Calculate stages
        stages_rect, stages_strip, optimal_feed, stage_data = calculate_stages(
            x_d, x_b, x_feed_intersect, y_feed_intersect, 
            slope_rect, intercept_rect, slope_strip, intercept_strip
        )
        
        total_stages = stages_rect + stages_strip
        
        # Flow rates
        V = D * (R + 1)  # Vapor rate in rectifying section
        L = R * D  # Liquid rate in rectifying section
        
        return {
            'F': F, 'D': D, 'B': B,
            'x_f': x_f, 'x_d': x_d, 'x_b': x_b,
            'q': q, 'R': R, 'R_min': R_min, 'R_factor': R_factor,
            'slope_rect': slope_rect, 'intercept_rect': intercept_rect,
            'slope_strip': slope_strip, 'intercept_strip': intercept_strip,
            'x_feed_intersect': x_feed_intersect, 'y_feed_intersect': y_feed_intersect,
            'x_intersect': x_intersect, 'y_intersect': y_intersect,
            'stages_rect': stages_rect, 'stages_strip': stages_strip,
            'total_stages': total_stages, 'optimal_feed': optimal_feed,
            'V': V, 'L': L, 'stage_data': stage_data
        }
    
    def calculate_stages(x_d, x_b, x_feed_intersect, y_feed_intersect, 
                        slope_rect, intercept_rect, slope_strip, intercept_strip):
        """
        Calculate number of stages using McCabe-Thiele graphical method
        """
        stages_data = []
        
        # Start from distillate composition
        x_current = x_d
        y_current = x_d  # At total condenser
        stage_num = 0
        stages_rect = 0
        stages_strip = 0
        feed_stage = None
        
        # Rectifying section
        while x_current > x_feed_intersect and stage_num < 50:
            stage_num += 1
            stages_rect += 1
            
            # Step down to equilibrium curve
            x_new = x_from_y_equilibrium(y_current)
            
            # Step across to operating line
            if x_new >= x_feed_intersect:
                y_new = slope_rect * x_new + intercept_rect
                section = "Rectifying"
            else:
                # Switch to stripping section
                y_new = slope_strip * x_new + intercept_strip
                section = "Stripping"
                feed_stage = stage_num
                stages_rect -= 1
                stages_strip += 1
            
            stages_data.append({
                'Stage': stage_num,
                'Section': section,
                'x_liquid': x_new,
                'y_vapor': y_current,
                'x_next': x_current,
                'y_next': y_new
            })
            
            x_current = x_new
            y_current = y_new
            
            if x_current <= x_b * 1.01:  # Close to bottoms composition
                break
        
        # Continue with stripping section if not already switched
        if feed_stage is None:
            feed_stage = stages_rect
        
        while x_current > x_b * 1.01 and stage_num < 50:
            stage_num += 1
            stages_strip += 1
            
            # Step down to equilibrium curve
            x_new = x_from_y_equilibrium(y_current)
            
            # Step across to stripping operating line
            y_new = slope_strip * x_new + intercept_strip
            
            stages_data.append({
                'Stage': stage_num,
                'Section': "Stripping",
                'x_liquid': x_new,
                'y_vapor': y_current,
                'x_next': x_current,
                'y_next': y_new
            })
            
            x_current = x_new
            y_current = y_new
            
            if x_current <= x_b * 1.01:
                break
        
        return stages_rect, stages_strip, feed_stage, stages_data
    
    @render.plot
    def mccabe_thiele_plot():
        design = calculate_design()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Equilibrium curve
        x_eq = np.linspace(0, 1, 100)
        y_eq = [equilibrium_curve(x) for x in x_eq]
        ax.plot(x_eq, y_eq, 'b-', linewidth=2, label='Equilibrium Curve')
        
        # 45-degree line
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='y = x')
        
        # Operating lines
        x_op = np.linspace(0, 1, 100)
        
        # Rectifying section
        y_rect = design['slope_rect'] * x_op + design['intercept_rect']
        mask_rect = (x_op >= design['x_feed_intersect']) & (x_op <= design['x_d'])
        ax.plot(x_op[mask_rect], y_rect[mask_rect], 'r-', linewidth=2, label='Rectifying Operating Line')
        
        # Stripping section
        y_strip = design['slope_strip'] * x_op + design['intercept_strip']
        mask_strip = (x_op >= design['x_b']) & (x_op <= design['x_feed_intersect'])
        ax.plot(x_op[mask_strip], y_strip[mask_strip], 'g-', linewidth=2, label='Stripping Operating Line')
        
        # Feed line
        if design['q'] != 1:
            y_feed = design['q']/(design['q']-1) * x_op - design['x_f']/(design['q']-1)
            mask_feed = (x_op >= min(design['x_b'], design['x_feed_intersect'])) & (x_op <= max(design['x_d'], design['x_feed_intersect']))
            ax.plot(x_op[mask_feed], y_feed[mask_feed], 'm--', linewidth=2, label='Feed Line')
        
        # Key points
        ax.plot(design['x_d'], design['x_d'], 'ro', markersize=8, label=f'Distillate ({design["x_d"]:.3f})')
        ax.plot(design['x_b'], design['x_b'], 'go', markersize=8, label=f'Bottoms ({design["x_b"]:.3f})')
        ax.plot(design['x_f'], design['x_f'], 'mo', markersize=8, label=f'Feed ({design["x_f"]:.3f})')
        ax.plot(design['x_feed_intersect'], design['y_feed_intersect'], 'ko', markersize=6, label='Feed Stage')
        
        # Draw stages
        stage_data = design['stage_data']
        for i, stage in enumerate(stage_data):
            if i < len(stage_data) - 1:  # Don't draw last incomplete stage
                # Vertical line (equilibrium step)
                ax.plot([stage['x_next'], stage['x_liquid']], [stage['y_vapor'], stage['y_vapor']], 
                       'k-', alpha=0.7, linewidth=1)
                # Horizontal line (operating line step)
                ax.plot([stage['x_liquid'], stage['x_liquid']], [stage['y_vapor'], stage['y_next']], 
                       'k-', alpha=0.7, linewidth=1)
        
        ax.set_xlabel('x (liquid mole fraction acetone)', fontsize=12)
        ax.set_ylabel('y (vapor mole fraction acetone)', fontsize=12)
        ax.set_title('McCabe-Thiele Diagram for Acetone-Ethanol Distillation', fontsize=14)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        return fig
    
    @render.text
    def design_results():
        design = calculate_design()
        
        results = f"""DESIGN RESULTS:
═══════════════

Number of Stages:
- Rectifying Section: {design['stages_rect']} stages
- Stripping Section: {design['stages_strip']} stages
- Total Theoretical Stages: {design['total_stages']} stages
- Optimal Feed Stage: {design['optimal_feed']}

Reflux Ratio:
- Minimum Reflux Ratio (Rmin): {design['R_min']:.3f}
- Operating Reflux Ratio (R): {design['R']:.3f}
- Reflux Ratio Factor: {design['R_factor']:.2f}

Operating Conditions:
- q-factor: {design['q']:.3f}
- Feed thermal condition: {'Subcooled liquid' if design['q'] > 1 else 'Saturated liquid' if design['q'] == 1 else 'Partially vaporized' if design['q'] > 0 else 'Saturated vapor' if design['q'] == 0 else 'Superheated vapor'}

Tray Efficiency:
- Assuming 100% efficiency
- Number of actual trays = {design['total_stages']} trays
"""
        return results
    
    @render.text
    def material_balance():
        design = calculate_design()
        
        balance = f"""MATERIAL BALANCE:
═══════════════

Overall Balance:
- Feed Rate (F): {design['F']:.1f} kmol/hr
- Distillate Rate (D): {design['D']:.1f} kmol/hr
- Bottoms Rate (B): {design['B']:.1f} kmol/hr

Component Balance (Acetone):
- Feed: {design['F']:.1f} × {design['x_f']:.3f} = {design['F'] * design['x_f']:.1f} kmol/hr
- Distillate: {design['D']:.1f} × {design['x_d']:.3f} = {design['D'] * design['x_d']:.1f} kmol/hr
- Bottoms: {design['B']:.1f} × {design['x_b']:.3f} = {design['B'] * design['x_b']:.1f} kmol/hr

Internal Flow Rates:
- Vapor rate in rectifying (V): {design['V']:.1f} kmol/hr
- Liquid rate in rectifying (L): {design['L']:.1f} kmol/hr
- Liquid rate in stripping: {design['L'] + design['F']:.1f} kmol/hr
- Vapor rate in stripping: {design['V']:.1f} kmol/hr

Recovery:
- Acetone recovery in distillate: {(design['D'] * design['x_d']) / (design['F'] * design['x_f']) * 100:.1f}%
- Ethanol recovery in bottoms: {(design['B'] * (1 - design['x_b'])) / (design['F'] * (1 - design['x_f'])) * 100:.1f}%
"""
        return balance
    
    @render.table
    def stage_analysis():
        design = calculate_design()
        stage_data = design['stage_data']
        
        if not stage_data:
            return pd.DataFrame()
        
        # Convert to DataFrame
        df = pd.DataFrame(stage_data)
        
        # Round numerical values
        numerical_cols = ['x_liquid', 'y_vapor', 'x_next', 'y_next']
        for col in numerical_cols:
            if col in df.columns:
                df[col] = df[col].round(4)
        
        # Rename columns for better display
        df = df.rename(columns={
            'Stage': 'Stage #',
            'Section': 'Section',
            'x_liquid': 'x (liquid)',
            'y_vapor': 'y (vapor)',
            'x_next': 'x (next)',
            'y_next': 'y (next)'
        })
        
        return df

app = App(app_ui, server)

if __name__ == "__main__":
    app.run(debug=True)