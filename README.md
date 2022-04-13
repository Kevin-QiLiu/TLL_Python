# README
This respoitory includes all python source code for Qi Liu's Master project in semiconductor physics group, Cavendish laboratory, University of Cambridge.

Contact Qi Liu at ql315@cam.ac.uk or qi.liu.apply@outlook.com 

## Wafer simulation
This file includes all the simulations done with the new wafer where width of the barrier is doubled but the tunnelling probability maintain the same as existing wafers. 
The Python module nexntnanopy is imported as an essential tool. 
For more information, see [nextnanopy's tutorial](https://github.com/nextnanopy/nextnanopy/blob/master/README.md):

## length variable device simulation
This file includes simulations for new device design, where we can vary the 1d channel length <em>in situ</em>.

## TODO
* There is some size problem in the method ```get_channel_length```
```python
    def get_channel_length(cls, matrix, edge_axis, slice_axis, tol, x_0, x_1):
        """return effective x upper limit for conducting channel, the channel length and width as a function of x"""
        global x_upper_limit
        num = 100
        width = []
        lim = np.linspace(x_0, x_1, num)
        count = 0
        for slice_pos in lim:
            channel_temp = cls(matrix, slice_pos, edge_axis, slice_axis)
            width_temp = channel_temp.get_channel_width()
            width.append(width_temp)
            if width_temp > tol:
                x_upper_limit = slice_pos
                break
            # elif width_temp == 1000:
            #     x_upper_limit = x_1
            count += 1
        return x_upper_limit, x_upper_limit - x_0, lim, width, count

```
## History of changes
### April 13th, 2022
* Initial upload



