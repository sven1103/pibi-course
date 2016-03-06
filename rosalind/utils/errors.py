class WrongAlphabetError(Exception):
    def __init__(self, *args):
        Exception.__init__(self, "Letter {0} is not in the {1} alphabet!"
                           .format(args[0], args[1]))


