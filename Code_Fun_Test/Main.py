from Function import *

my_car = car(False, 0, "Toyota")

my_car.increase_speed(20)  # This will print a message indicating that the car is not started.

baby_car = my_car.have_baby("Baby Honda")  # This will create a new baby car with the model "Honda" and print a message.

print(baby_car.get_info())  # This will return the information about the baby car's model.

print(my_car.get_info())  # This will return the information about the baby car's model.

baby_car.start()  # This will start the baby car and print a message.
baby_car.increase_speed(30)  # This will increase the speed of the baby car by 30 km/h and print the new speed.