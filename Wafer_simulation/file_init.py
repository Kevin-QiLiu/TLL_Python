import nextnanopy as nn

class INPUT_FILE:

    def __init__(self, path):
        self.path = path
        self.epath = str(path.format('electron_density'))