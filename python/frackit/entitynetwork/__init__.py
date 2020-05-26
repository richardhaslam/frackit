from ._entitynetwork import *

# python implementation of the constraints class that accepts evaluation on lists
class EntityNetworkConstraints(_entitynetwork._EntityNetworkConstraints):

    def evaluate(self, a, b):
        """
        Evaluates the constraints between the arguments a and b, both being
        either an instance of a geometry class or a list of geometry classes.

        Parameters:
        a: a geometry or list of geometries
        b: a geometry or list of geometries
        """
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


# python implementation of the entity network builder class
class EntityNetworkBuilder(_entitynetwork._EntityNetworkBuilderWrapper):

    def addEntities(self, entities):
        """
        adds all the entities in the provided list to the network.
        """
        for entity in entities: super().addEntity(entity)

    def addSubDomainEntities(self, entities, subDomainId):
        """
        adds all the given entities to be embedded in the subdomain with the given id.
        """
        for entity in entities: super().addSubDomainEntity(entity, subDomainId)

# python implementation of the contained entity network builder class
class ContainedEntityNetworkBuilder(_entitynetwork._ContainedEntityNetworkBuilderWrapper):

    def addSubDomainEntities(self, entities, subDomainId):
        """
        adds all the given entities to be embedded in the subdomain with the given id.
        """
        for entity in entities: super().addSubDomainEntity(entity, subDomainId)
