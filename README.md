# Power-system-state-estimation-using-the-gradient-descent-optimization
## Power System State Estimation – 2-Bus Model

This project implements a simplified **state estimation** method for a power system with **two buses**. The goal is to estimate the system state using real and reactive power flow measurements.

### System Overview
The **system state** consists of three unknowns:
- `V1`: Voltage magnitude at Bus 1  
- `V2`: Voltage magnitude at Bus 2  
- `θ2`: Voltage phase angle at Bus 2 (Bus 1 angle `θ1` is fixed at 0 as reference)

### 📐 Measurement Equations
Power flow measurements in both directions between the two buses:
- `P12`, `Q12`: Real and reactive power from Bus 1 to Bus 2  
- `P21`, `Q21`: Real and reactive power from Bus 2 to Bus 1  

These are modeled as nonlinear functions of the state variables:
