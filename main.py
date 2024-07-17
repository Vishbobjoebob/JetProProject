from designer import engineDesigner

engine = engineDesigner(220, 11.0, 1.10, 15, 1.2, 2.0, 0.06, 0.021, 0.0040, 0.0050)

engine.diffuser()
engine.bypass_fan()
engine.compressor()
engine.combustor_main()
engine.turbine()
engine.mixer_turbine()
engine.combustor_interturbine()
engine.turbine_fan()
engine.combustor_afterburner()
engine.nozzle_core()
engine.nozzle_fan()
engine.mixer_nozzle()
engine.nozzle_combined()

output = [engine.T04, engine.p04, engine.Tmax_t]

print(output)