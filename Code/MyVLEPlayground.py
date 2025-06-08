from shiny import App, render, ui, reactive
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import pandas as pd
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import io
import base64

# Define Antoine equation parameters for n-hexane and n-heptane
def psat1(temp):
    """Vapor pressure of n-hexane (bar)"""
    return 10**(4.00266 - 1171.53/(temp + 224.216))

def psat2(temp):
    """Vapor pressure of n-heptane (bar)"""
    return 10**(4.04867 - 1355.126/(temp + 209.367))

def px_function(x, temp):
    """Bubble point pressure"""
    return x * psat1(temp) + (1 - x) * psat2(temp)

def py_function(x, temp):
    """Dew point pressure"""
    return 1 / (x/psat1(temp) + (1-x)/psat2(temp))

def generate_equilibrium_data():
    """Generate equilibrium data for T-x-y diagram"""
    x_vals = np.linspace(0, 1, 101)
    t_bubble = []
    t_dew = []
    y_vals = []
    
    for x in x_vals:
        if x == 0:
            # Pure component 2
            t_bubble.append(98.42)  # Boiling point of n-heptane at 1 bar
            t_dew.append(98.42)
            y_vals.append(0)
        elif x == 1:
            # Pure component 1
            t_bubble.append(68.73)  # Boiling point of n-hexane at 1 bar
            t_dew.append(68.73)
            y_vals.append(1)
        else:
            # Find bubble point temperature
            def bubble_eq(T):
                return px_function(x, T) - 1.0  # 1 bar
            
            try:
                T_bubble = fsolve(bubble_eq, 80)[0]
                t_bubble.append(T_bubble)
                
                # Calculate y from bubble point
                y = x * psat1(T_bubble) / 1.0
                y_vals.append(y)
                
                # Find dew point temperature
                def dew_eq(T):
                    return py_function(x, T) - 1.0  # 1 bar
                
                T_dew = fsolve(dew_eq, 80)[0]
                t_dew.append(T_dew)
                
            except:
                t_bubble.append(80)
                t_dew.append(80)
                y_vals.append(x)
    
    return x_vals, t_bubble, t_dew, y_vals

# Generate equilibrium data
x_eq, t_bubble, t_dew, y_eq = generate_equilibrium_data()

# Create interpolation functions
tx_interp = interp1d(x_eq, t_bubble, kind='linear', fill_value='extrapolate')
ty_interp = interp1d(y_eq, t_dew, kind='linear', fill_value='extrapolate')

app_ui = ui.page_fluid(
    ui.h2("Binary Mixture VLE Interactive Tool"),
    ui.br(),
    
    ui.row(
        ui.column(3,
            ui.input_radio_buttons(
                "diagram_type", 
                "Diagram Type:",
                choices={"pxy": "P-x-y", "txy": "T-x-y"},
                selected="pxy"
            ),
            ui.br(),
            ui.panel_conditional(
                "input.diagram_type === 'pxy'",
                ui.input_slider(
                    "temperature", 
                    "Temperature (°C):", 
                    min=85, max=125, value=115, step=1
                ),
                ui.input_slider(
                    "pressure_point", 
                    "Pressure (bar):", 
                    min=0.5, max=4.0, value=1.5, step=0.1
                ),
                ui.input_slider(
                    "x_composition", 
                    "n-hexane mole fraction:", 
                    min=0.01, max=0.99, value=0.45, step=0.01
                )
            ),
            ui.panel_conditional(
                "input.diagram_type === 'txy'",
                ui.input_slider(
                    "pressure", 
                    "Pressure (bar):", 
                    min=0.5, max=2.0, value=1.0, step=0.1
                ),
                ui.input_slider(
                    "temperature_point", 
                    "Temperature (°C):", 
                    min=45, max=155, value=115, step=1
                ),
                ui.input_slider(
                    "x_composition_txy", 
                    "n-hexane mole fraction:", 
                    min=0.01, max=0.99, value=0.45, step=0.01
                )
            )
        ),
        
        ui.column(9,
            ui.row(
                ui.column(8,
                    ui.output_plot("phase_diagram", height="500px")
                ),
                ui.column(4,
                    ui.output_plot("lever_rule", height="500px")
                )
            )
        )
    ),
    
    ui.br(),
    ui.row(
        ui.column(12,
            ui.output_text_verbatim("phase_info")
        )
    )
)

def server(input, output, session):
    
    @reactive.calc
    def calculate_phase_equilibrium():
        if input.diagram_type() == "pxy":
            T = input.temperature()
            P_point = input.pressure_point()
            z = input.x_composition()
            
            # Calculate bubble and dew point pressures
            P_bubble = px_function(z, T)
            P_dew = py_function(z, T)
            
            # Determine phase and compositions
            if P_point < P_dew:
                # Vapor phase only
                phase = "vapor"
                x_liquid = z
                y_vapor = z
                vapor_fraction = 1.0
            elif P_point > P_bubble:
                # Liquid phase only
                phase = "liquid"
                x_liquid = z
                y_vapor = z
                vapor_fraction = 0.0
            else:
                # Two-phase region
                phase = "two-phase"
                # Solve for liquid composition
                def liquid_eq(x):
                    return px_function(x, T) - P_point
                x_liquid = fsolve(liquid_eq, z)[0]
                
                # Solve for vapor composition
                def vapor_eq(x):
                    return py_function(x, T) - P_point
                y_vapor = fsolve(vapor_eq, z)[0]
                
                # Calculate vapor fraction using lever rule
                vapor_fraction = (y_vapor - z) / (y_vapor - x_liquid)
            
            return {
                'phase': phase,
                'x_liquid': x_liquid,
                'y_vapor': y_vapor,
                'vapor_fraction': vapor_fraction,
                'P_bubble': P_bubble,
                'P_dew': P_dew,
                'T': T,
                'P_point': P_point,
                'z': z
            }
            
        else:  # txy diagram
            P = input.pressure()
            T_point = input.temperature_point()
            z = input.x_composition_txy()
            
            # For T-x-y at constant pressure, use equilibrium data
            T_bubble_z = tx_interp(z)
            T_dew_z = ty_interp(z) if z < 1 else tx_interp(z)
            
            # Determine phase and compositions
            if T_point > T_dew_z:
                # Vapor phase only
                phase = "vapor"
                x_liquid = z
                y_vapor = z
                vapor_fraction = 1.0
            elif T_point < T_bubble_z:
                # Liquid phase only
                phase = "liquid"
                x_liquid = z
                y_vapor = z
                vapor_fraction = 0.0
            else:
                # Two-phase region
                phase = "two-phase"
                # Interpolate to find compositions
                # Find x where T_bubble(x) = T_point
                x_range = np.linspace(0, 1, 1000)
                t_range = [tx_interp(x) for x in x_range]
                x_liquid = x_range[np.argmin(np.abs(np.array(t_range) - T_point))]
                
                # Find y where T_dew(y) = T_point
                y_range = np.linspace(0, 1, 1000)
                t_dew_range = [ty_interp(y) for y in y_range]
                y_vapor = y_range[np.argmin(np.abs(np.array(t_dew_range) - T_point))]
                
                # Calculate vapor fraction using lever rule
                if y_vapor != x_liquid:
                    vapor_fraction = (y_vapor - z) / (y_vapor - x_liquid)
                else:
                    vapor_fraction = 0.5
            
            return {
                'phase': phase,
                'x_liquid': x_liquid,
                'y_vapor': y_vapor,
                'vapor_fraction': vapor_fraction,
                'T_bubble': T_bubble_z,
                'T_dew': T_dew_z,
                'T_point': T_point,
                'P': P,
                'z': z
            }
    
    @render.plot
    def phase_diagram():
        equilibrium = calculate_phase_equilibrium()
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        if input.diagram_type() == "pxy":
            T = equilibrium['T']
            x_range = np.linspace(0, 1, 100)
            p_bubble = [px_function(x, T) for x in x_range]
            p_dew = [py_function(x, T) for x in x_range]
            
            ax.plot(x_range, p_bubble, 'b-', linewidth=2, label='Bubble point')
            ax.plot(x_range, p_dew, 'g-', linewidth=2, label='Dew point')
            
            # Plot operating point and tie line
            z = equilibrium['z']
            P_point = equilibrium['P_point']
            
            if equilibrium['phase'] == "two-phase":
                x_liq = equilibrium['x_liquid']
                y_vap = equilibrium['y_vapor']
                
                # Tie line
                ax.plot([x_liq, y_vap], [P_point, P_point], 'k--', linewidth=2)
                ax.plot([x_liq, x_liq], [0, P_point], 'b:', linewidth=2)
                ax.plot([y_vap, y_vap], [0, P_point], 'g:', linewidth=2)
                ax.plot([z, z], [0, P_point], 'r:', linewidth=2)
            else:
                ax.plot([z, z], [0, P_point], 'r:', linewidth=2)
            
            ax.plot(z, P_point, 'ro', markersize=8, label='Operating point')
            
            ax.set_xlabel('n-hexane mole fraction')
            ax.set_ylabel('Pressure (bar)')
            ax.set_title(f'P-x-y diagram at {T}°C')
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 4.5)
            ax.text(0.15, 3.8, 'liquid', fontsize=14, color='gray')
            ax.text(0.85, 0.5, 'vapor', fontsize=14, color='gray')
            
        else:  # txy diagram
            ax.plot(x_eq, t_bubble, 'b-', linewidth=2, label='Bubble point')
            ax.plot(y_eq, t_dew, 'g-', linewidth=2, label='Dew point')
            
            # Plot operating point and tie line
            z = equilibrium['z']
            T_point = equilibrium['T_point']
            
            if equilibrium['phase'] == "two-phase":
                x_liq = equilibrium['x_liquid']
                y_vap = equilibrium['y_vapor']
                
                # Tie line
                ax.plot([x_liq, y_vap], [T_point, T_point], 'k--', linewidth=2)
                ax.plot([x_liq, x_liq], [40, T_point], 'b:', linewidth=2)
                ax.plot([y_vap, y_vap], [40, T_point], 'g:', linewidth=2)
                ax.plot([z, z], [40, T_point], 'r:', linewidth=2)
            else:
                ax.plot([z, z], [40, T_point], 'r:', linewidth=2)
            
            ax.plot(z, T_point, 'ro', markersize=8, label='Operating point')
            
            ax.set_xlabel('n-hexane mole fraction')
            ax.set_ylabel('Temperature (°C)')
            ax.set_title(f'T-x-y diagram at {equilibrium["P"]} bar')
            ax.set_xlim(0, 1)
            ax.set_ylim(40, 160)
            ax.text(0.15, 50, 'liquid', fontsize=14, color='gray')
            ax.text(0.85, 150, 'vapor', fontsize=14, color='gray')
        
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        return fig
    
    @render.plot
    def lever_rule():
        equilibrium = calculate_phase_equilibrium()
        
        fig, ax = plt.subplots(figsize=(3, 6))
        
        vapor_frac = equilibrium['vapor_fraction']
        liquid_frac = 1 - vapor_frac
        
        # Draw rectangles for liquid and vapor fractions
        if liquid_frac > 0:
            rect_liquid = Rectangle((0, 0), 0.5, liquid_frac, 
                                  facecolor='blue', alpha=0.7)
            ax.add_patch(rect_liquid)
            
            # Add text for liquid composition
            ax.text(0.25, liquid_frac/2, f'x_H = {equilibrium["x_liquid"]:.2f}', 
                   ha='center', va='center', color='white', fontsize=12, fontweight='bold')
        
        if vapor_frac > 0:
            rect_vapor = Rectangle((0, liquid_frac), 0.5, vapor_frac, 
                                 facecolor='green', alpha=0.7)
            ax.add_patch(rect_vapor)
            
            # Add text for vapor composition
            ax.text(0.25, liquid_frac + vapor_frac/2, f'y_H = {equilibrium["y_vapor"]:.2f}', 
                   ha='center', va='center', color='white', fontsize=12, fontweight='bold')
        
        ax.set_xlim(-0.1, 0.6)
        ax.set_ylim(0, 1)
        ax.set_ylabel('liquid and vapor amounts (mol)', rotation=90)
        ax.set_title('Phase Amounts')
        ax.set_xticks([])
        
        # Add scale on the left
        for i in np.arange(0, 1.1, 0.2):
            ax.axhline(y=i, color='gray', linestyle='-', alpha=0.3)
            ax.text(-0.05, i, f'{i:.1f}', ha='right', va='center')
        
        plt.tight_layout()
        return fig
    
    @render.text
    def phase_info():
        equilibrium = calculate_phase_equilibrium()
        
        info = f"Phase: {equilibrium['phase'].upper()}\n"
        info += f"Overall composition (z): {equilibrium['z']:.3f}\n"
        
        if equilibrium['phase'] == "two-phase":
            info += f"Liquid composition (x): {equilibrium['x_liquid']:.3f}\n"
            info += f"Vapor composition (y): {equilibrium['y_vapor']:.3f}\n"
            info += f"Vapor fraction: {equilibrium['vapor_fraction']:.3f}\n"
            info += f"Liquid fraction: {1-equilibrium['vapor_fraction']:.3f}\n"
        
        if input.diagram_type() == "pxy":
            info += f"\nBubble point pressure: {equilibrium['P_bubble']:.3f} bar\n"
            info += f"Dew point pressure: {equilibrium['P_dew']:.3f} bar\n"
        else:
            info += f"\nBubble point temperature: {equilibrium['T_bubble']:.1f} °C\n"
            info += f"Dew point temperature: {equilibrium['T_dew']:.1f} °C\n"
        
        return info

app = App(app_ui, server)

if __name__ == "__main__":
    app.run(debug=True)