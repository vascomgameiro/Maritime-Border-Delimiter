from daup2425p1 import Point, ValidPoints, Delimitation, FileHandler
class ConvexHull:
    """
    A class to calculate the convex hull from a set of valid points.
    
    The convex hull is the smallest convex shape that encloses a set of points.
    This class uses the orientation method to find the hull by iterating over the points, 
    ensuring that only those points that form part of the hull are retained.
    
    Attributes:
        __points (ValidPoints): The set of valid points from which the convex hull is calculated.
    
    Methods:
    find_delimitation():Returns the convex hull (delimitation) for the given set of points. 
    """

    def __init__(self, valid_points: ValidPoints):
        """
        Initializes the ConvexHull class with a set of valid points.
        
        Args:
            valid_points (ValidPoints): An instance containing the points from which the convex hull is to be computed.
        
        Complexity: O(1)
        """
        self.__points = valid_points

    @staticmethod
    def __orientation(p1: Point, p2: Point, p3: Point) -> int:
        """
        Determines the orientation of the triplet (p1, p2, p3).
        
        Given three points, the method calculates the orientation to decide if the triplet forms a clockwise, 
        counterclockwise, or collinear angle. This helps in determining whether to add or discard a point when 
        forming the convex hull.
        
        Args:
            p1 (Point): The first point of the triplet (previous point in the sequence).
            p2 (Point): The second point of the triplet (current point being evaluated).
            p3 (Point): The third point of the triplet (next point in the sequence).
        
        Returns:
            int: 
                - 0 if the points are collinear,
                - 1 if they are in clockwise orientation,
                - 2 if they are in counterclockwise orientation.
        
        Complexity: O(1)
        """
        val = (p2.get_longitude() - p1.get_longitude()) * (p3.get_latitude() - p2.get_latitude()) - \
              (p2.get_latitude() - p1.get_latitude()) * (p3.get_longitude() - p2.get_longitude())
        
        if val == 0:
            return 0  # collinear
        elif val > 0:
            return 1  # clockwise
        else:
            return 2  # counterclockwise

    def find_delimitation(self) -> Delimitation:
        """
        Finds and returns the convex hull (delimitation) for the given set of points.
        
        This method initializes the delimitation process by finding the starting point with the lowest y-coordinate, 
        and then iterates over the remaining points, using the orientation method to decide whether to add or remove points.
        
        Points are added or removed based on their angular orientation relative to the last two points in the hull.
        
        Returns:
            Delimitation: The object containing the points that make up the convex hull.
        
        Complexity: O(plog(p)), p being the number of points in the Delimitation object
        """
        d = Delimitation()
        points = self.__points.get_all_points()
        dictionary1 = dict()

        start = min(points, key=lambda p: p.get_longitude())

        for i in range(self.__points.get_size()):
            ponto_ref = Point("ref", start.get_latitude() - 1, start.get_longitude() - 1)
            angle = start.get_forward_angle(ponto_ref, points[i])
            dictionary1[points[i]] = angle

        dictionary1.pop(start)
        sorted_points = sorted(dictionary1, key=dictionary1.get)

        d.add_point(start)

        while len(sorted_points) > 0:
            if d.size() < 2:
                d.add_point(sorted_points[0])
                sorted_points.pop(0)
                continue

            (p2, p3) = d.get_last_two()
            
            curr_point = sorted_points[0]
            orientation_result = self.__orientation(p2, p3, curr_point)
            
            if orientation_result == 1:  
                d.pop_point()
            else:
                d.add_point(curr_point)  
                sorted_points.pop(0)
              
        return d
