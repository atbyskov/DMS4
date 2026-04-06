class bike:
    def __init__(self):
        self.start_Status = False
        self.speed = 0
        self.model = "Yamaha"
    
    def start(self):
        self.start_Status = True
        print("The bike has started.")
 
    def increase_speed(self, increment):
        if not self.start_Status:
            print("The bike is not started. Please start the bike before increasing speed.")
            return
        else:
            self.speed += increment
            print(f"The bike's speed has increased to {self.speed} km/h.")
    
    def make_baby(self, baby_model):
        baby_bike = bike()
        baby_bike.model = baby_model
        print(f"A new baby bike has been created with model {baby_model}.")
        return baby_bike

 # text = "Hej med dig"  string
 # ny_text = text.uppercase()