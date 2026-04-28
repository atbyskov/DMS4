import h5py
import numpy as np

with h5py.File('Optimization_Logs/pyslsqp_history.hdf5', 'r') as f:
    print("HDF5 Structure:")
    print(list(f.keys()))
    print()
    
    # Find the last iteration
    iter_keys = sorted([k for k in f.keys() if k.startswith('iter_')])
    if iter_keys:
        last_iter = iter_keys[-1]
        print(f"Last iteration: {last_iter}")
        print(f"Keys in {last_iter}: {list(f[last_iter].keys())}")
        print()
        
        # Get the best results from the last iteration
        last_x = f[last_iter]['x'][()]  # Best design variables
        last_obj = f[last_iter]['objective'][()]  # Best mass
        
        print("="*60)
        print("BEST OPTIMIZATION RESULT")
        print("="*60)
        print(f"\nOptimal Mass: {last_obj:.3f} kg")
        print(f"\nOptimal Design Variables:")
        
        # Variable names in order (from Main.py)
        var_names = ['rad', 'd0_1', 'd0_2', 'd0_3', 'd0_4', 'd0_5',
                     't0_1', 't0_2', 't0_3', 't0_4', 't0_5',
                     'd1_1', 'd1_2', 'd1_3', 'd1_4', 'd1_5',
                     't1_1', 't1_2', 't1_3', 't1_4', 't1_5']
        
        for name, value in zip(var_names, last_x):
            print(f"  {name:8s}: {value:10.4f}")
        print("="*60)
    else:
        print("No iterations found in file")