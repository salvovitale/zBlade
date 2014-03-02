from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Display.SimpleGui import init_display
 
display, start_display, add_menu, add_function_to_menu = init_display()
 
my_box = BRepPrimAPI_MakeBox(10.,20.,30.).Shape()
 
display.DisplayShape(my_box, update=True )
start_display()
#this is a stupid try out
