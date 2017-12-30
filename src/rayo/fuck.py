import matplotlib.pyplot as plt

# a 4x4 box (counterclockwise)
ext_x = [2, -2, -2, 2, 2]
ext_y = [2, 2, -2, -2, 2]

# a 2x2 hole in the box (clockwise)
int_x = [item/2.0 for item in ext_x][::-1]
int_y = [item/2.0 for item in ext_y][::-1]

# if you don't specify a color, you will see a seam
plt.fill(ext_x+int_x, ext_y+int_y, color='blue')

plt.show()
