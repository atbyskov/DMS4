class car:
    def __init__(self, start_Status: bool, speed: int, model: str):
        self.start_Status = start_Status
        self.speed = speed
        self.model = model

    def start(self):
        self.start_Status = True
        print("The car has started.")

    def increase_speed(self, increment):
        if not self.start_Status:
            print("The car is not started. Please start the car before increasing speed.")
            return
        else:
            self.speed += increment
            print(f"The car's speed has increased to {self.speed} km/h.")
            # self.speed = self.speed + increment

    def have_baby(self, baby_model):
        baby_car = car(False, 0, baby_model)
        print(f"A new baby car has been created with model {baby_model}.")
        return baby_car

    def stop(self):
        self.speed = 0
        self.start_Status = False
        print("The car has stopped.")

    def get_info(self):
        return f"The model of the car is a {self.model}"
    
    def change_model(self, new_model):
        self.model = new_model
        print(f"The car model has been changed to {self.model}.")