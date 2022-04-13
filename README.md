# README
This respoitory includes all python source code for Qi Liu's Master project in Semiconductor Physics group, Cavendish Laboratory, University of Cambridge.

Contact Qi Liu at ql315@cam.ac.uk or qi.liu.apply@outlook.com 

## Wafer simulation
This file includes all the simulations done with the new wafer where width of the barrier is doubled but the tunnelling probability maintain the same as existing wafers. 
The Python module nexntnanopy is imported as an essential tool. 
For more information, see [nextnanopy's tutorial](https://github.com/nextnanopy/nextnanopy/blob/master/README.md):

## length variable device simulation
This file includes simulations for new device design, where we can vary the 1d channel length <em>in situ</em>.

## TODO
* Optimise grid points: reduce grid point at central of channel and also regrowth interface
* Try single well model where tunnelling probability is same to that as in superlattice
* Adjust gate geometry to compensate potential tail
* Non-equibilirium Green's function package: ```NextNano.NEGF```
## History of changes
### April 13th, 2022
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



