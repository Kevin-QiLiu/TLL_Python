# README
This respoitory includes all python source code for Qi Liu's Master project in Semiconductor Physics group, Cavendish Laboratory, University of Cambridge.

Contact Qi Liu at ql315@cam.ac.uk or qi.liu.apply@outlook.com 

## Wafer simulation
This file includes all the simulations done with the new wafer where width of the barrier is doubled but the tunnelling probability maintain the same as existing wafers. 
The Python module nexntnanopy is imported as an essential tool. 
For more information, see [nextnanopy's tutorial](https://github.com/nextnanopy/nextnanopy/blob/master/README.md):

## Length-variable device simulation
This file includes simulations for new device design, where we can vary the 1d channel length <em>in situ</em>.

## TODO
* Vary gate voltages using Davies's method to test channel length, develop a 
* Optimise nextnano grid points: reduce grid point at central of channel and also regrowth interface (pirority)
* Try single well model where tunnelling probability is same to that as in superlattice (less sugessted)
* Adjust gate geometry to compensate potential tail
* Try Non-equibilirium Green's function package: ```NextNano.NEGF```
## History of changes
### April 13th, 2022
* introduce a ```delta_tolerance``` so that the cursor can scan through the whole region without interupting 
```python
    @classmethod
    def get_channel_length(cls, matrix, x_axis, y_axis, tol, tol_delta, x_0, x_1):
        """return effective x upper/ lower limit and length for conducting channel,
        assuming that we have tolerance width = fac * width at midpoint to form 1d channel, the neighbourhood of
        tolerance is tol_delta, choice of it depends on grid point resolution"""
        itr_num = 2000
        endpts_set = []
        lim = np.linspace(x_0, x_1, itr_num)
        Qdot_check = cls(matrix, 0, y_axis, x_axis)
        if len(Qdot_check.get_Ef_intersect()[4]) >= 3:
            print('In Quantum Dot Regime! Check y_sliced potential! ')
        for slice_pos in lim:
            channel_temp = cls(matrix, slice_pos, x_axis, y_axis)
            if len(channel_temp.get_Ef_intersect()[2]) != 2:
                continue
            else:
                width_temp = channel_temp.get_channel_width()
                if abs(width_temp - tol) < tol_delta:
                    endpts_set.append(slice_pos)
        if len(endpts_set) == 0:
            raise Exception('No end point is chosen, increase tol_delta!')
        channel_x_1 = max(endpts_set)
        channel_x_0 = min(endpts_set)
        channel_length = channel_x_1 - channel_x_0
        return channel_length, channel_x_0, channel_x_1, endpts_set

```
* separate ```get_channel_length``` and ```get_vary_width``` method to resolve the size boundary issue.
```python
    @classmethod
    def get_vary_width(cls, matrix, edge_axis, slice_axis, x_0, x_1):
        width = []
        num = 100
        lim = np.linspace(x_0, x_1, num)
        for slice_pos in lim:
            channel_temp = cls(matrix, slice_pos, edge_axis, slice_axis)
            width_temp = channel_temp.get_channel_width()
            width.append(width_temp)
        return lim, np.asarray(width)
```
* Initial upload




