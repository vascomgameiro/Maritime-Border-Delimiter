from pythonds3 import Graph, Vertex

from daup2425p2 import Delimitation, ValidPoints


class OptimalDelimitation:
    """
    The OptimalDelimitation class is responsible for finding the delimitation with the maximum area
    from a set of valid points, given a distance constraint. It builds a graph of points and uses 
    Depth-First Search (DFS) to explore all possible delimitations (polygonal shapes), storing valid 
    delimitations by their calculated area. The delimitation with the maximum area is returned.

    Attributes:__distance : int,  __valid_points : ValidPoints

    Methods: find_delimitation(): Returns the optimal delimitation (Delimitation) for the given set of points. 
    """
    def __init__(self, valid_points: ValidPoints, distance: int):
        """
        This method initializes an OptimalDelimitation object.

        Arguments: valid_points(ValidPoints), distance (int)

        Complexity: O(1)
        """
        self.__distance = distance
        self.__valid_points = valid_points

    
    def __dfs_visit(self, current_vertex: Vertex, delimitation: Delimitation, area_registry: dict):
        """
        Performs a Depth-First Search (DFS) traversal starting from the given vertex, building
        possible delimitations and calculating their areas. The valid delimitations are stored
        in area_registry.

        Arguments: current_vertex: The vertex from which DFS is initiated (a Vertex object).
                   delimitation: A Delimitation object that tracks the sequence of points in the current path.
                   area_registry: A Dictionary to store valid delimitation areas mapped to the corresponding Delimitation object.
        Complexity:
        """
        current_vertex.color = "grey"
        for next_vertex in current_vertex.get_neighbors():
            if next_vertex.color == "white":
                delimitation.add_point(next_vertex.get_key())  
                if delimitation.size() >= 3: 
                    try:
                        total_area = delimitation.get_area()
                        area_registry[total_area] = delimitation.copy()
                    except ValueError: 
                        delimitation.pop_point()
                        continue
                self.__dfs_visit(next_vertex, delimitation, area_registry)
                if delimitation.size() > 0:
                    delimitation.pop_point()
        current_vertex.color = "white"


    
    def find_delimitation(self) -> Delimitation:
        """
        Finds and returns the delimitation with the maximum area by building a graph of points and performing DFS.

        The method constructs a graph where points are vertices and edges represent the valid connections
        between points based on the distance threshold. DFS is performed on each point, and the valid 
        delimitations are stored. The delimitation with the maximum area is returned.

        Returns: Delimitation object

        Complexity: O(p^2), p being the number of points in the ValidPoints object

        """
        area_registry = dict()
        points = self.__valid_points.get_all_points()        
        g = Graph()
        for i in range(self.__valid_points.get_size()):
            neigbors = self.__valid_points.get_points_vicinity(points[i], self.__distance)
            for neighbor in neigbors:
                if points[i] != neighbor:
                    g.add_edge(points[i], neighbor, points[i].distance(neighbor))
        
        for vertex in g.get_vertices():
            if g.get_vertex(vertex).color != "black":
                g.get_vertex(vertex).color = "white"
            d = Delimitation()
            d.add_point(vertex)  
            self.__dfs_visit(g.get_vertex(vertex), d, area_registry)
            g.get_vertex(vertex).color = "black"
        
        valid_area_registry = {area: delimitation for area, delimitation in area_registry.items() if delimitation.is_valid_delimitation(self.__valid_points, self.__distance)}
        if len(valid_area_registry) > 0:
            maximum = max(valid_area_registry)
            return area_registry[maximum]
        return Delimitation()
