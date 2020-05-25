from ._entitynetwork import *

class EntityNetworkConstraints(_entitynetwork._EntityNetworkConstraints):

    def evaluate(self, a, b):

        isListA = isinstance(a, list)
        isListB = isinstance(b, list)

        def evalList(theList, geo):
            return all( super(EntityNetworkConstraints, self).evaluate(item, b) for item in theList )

        if isListA and not isListB:
            return evalList(a, b)
        elif isListB and not isListA:
            return evalList(b, a)
        elif isListA and isListB:
            return all( evalList(b, item) for item in a )
        # here we assume the two arguments to be geometries
        else:
            return super().evaluate(a, b)
