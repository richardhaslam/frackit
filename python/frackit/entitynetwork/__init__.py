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

# python implementation of the constraints matrix class
class EntityNetworkConstraintsMatrix:

    def __init__(self):
        self.constraintsMatrix = {}

    def addConstraints(self, constraints, idPairs):

        def addIdPair(pair):
            if len(pair) != 2:
                raise RuntimeError("Id pairs are expected to be lists of length 2")

            id1 = pair[0]
            id2 = pair[1]
            if self.constraintsMatrix.get(id1) == None:
                self.constraintsMatrix[id1] = {id2 : [constraints]}
            elif self.constraintsMatrix.get(id1).get(id2) == None:
                self.constraintsMatrix[id1][id2] = [constraints]
            else:
                raise RuntimeError("Id pair already taken!")

        if not isinstance(idPairs, list):
            raise RuntimeError("This expects id pairs in the form of index lists")
        if not isinstance(idPairs[0], list): # single id pair
            addIdPair(idPairs)
        else:
            for pair in idPairs:
                addIdPair(pair)

    def evaluate(self, entitySets, entity, id):
        """
        Evaluates the constraints for all entity sets that have constraints defined
        against the entity set with the given id.

        Parameters:
        entitySets: dictionary that maps ids to lists of entities
        entity: the entity candidate for which the constraints are to be checked
        id: the id of the set to which the entity candidate belongs.
        """
        for setId in self.constraintsMatrix:

            entitySet = entitySets.get(setId)

            # set is empty, skip rest
            if entitySet == None:
                continue

            elif self.constraintsMatrix.get(setId).get(id) != None:
                # a constraint is defined between this set and the given id
                # evaluate all constraints defined for this pair
                constraints = self.constraintsMatrix.get(setId).get(id)
                if any(not c.evaluate(entitySet, entity) for c in constraints):
                    return False
        return True



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
