import matplotlib.pyplot as plt
from src.point import Point
from src.valid_points import ValidPoint

class Delimitation:
    """
    The Delimitation ADT models a polygon that can be constructed by sequentially adding 
    points.
 
    Attributes:
        __points (list): A list that stores all the points in the order they are added.
        __first_point (Point or None): The first point added to the Delimitation.
        __third_last_point (Point or None): The third-to-last point added to the Delimitation.
        __second_last_point (Point or None): The second-to-last point added to the Delimitation.
        __last_point (Point or None): The last point added to the Delimitation.
    
    Methods:
        add_point(Point): Adds a point to the Delimitation and updates relevant attributes.
        get_first(): Returns the first point added to the Delimitation, or None if no points exist.
        get_last_two(): Returns the last two points added as a tuple.
        get_points(): Returns the list of all points in the order they were added.
        pop_point(): Removes and returns the last point added to the Delimitation.
        get_area(): Calculates and returns the area of the polygon.
        copy(): Returns a deep copy of the current Delimitation object.
        intersects(Point, Point, Point, Point): Determines whether two line segments formed by the given points intersect.
        crosses_delimitation(Point, Point): Checks if a line segment formed by two points intersects with any existing segments in the Delimitation.
        show(ValidPoint): Displays the current Delimitation and points using a plotting tool.
    """
    def __init__(self):
        """
        This method initializes an empty Delimitation object.

        Complexity: O(1)
        """
        self.__points=[]
        self.__first_point=None
        self.__last_point=None
        self.__second_last_point=None
        self.__third_last_point=None
    
    def get_points(self) -> list:
        """
        This method retrieves the Points that belong to the Delimitation object by the order
        they were added.

        Returns: list of Point objects

        Complexity: O(1)
        """
        return self.__points
     
    def get_first(self) -> Point:
        """
        This method retrieves the first point that was added to the Delimitation object.

        Returns: Point object

        Complexity: O(1)
        """
        return self.__first_point
    
    def get_last_two(self) -> tuple:
        """ 
        This method retrieves the last two Point objects that have been added to the
        Delimitation object,the rightmost being the last one added.
        
        Returns: tuple(Point, Point)

        Complexity: O(1)
        """
        return (self.__second_last_point, self.__last_point)
    
    def get_area(self) -> float:
        """
        This method calculates the area of the polygon that is represented. If
        the Delimitation object has less than 3 points, the method raises a ValueError.
        
        Returns: float

        Complexity: O(p), p being the number of points
        """
        vertices=self.get_points()
        n = len(vertices)
        # if n < 3:
        #     raise ValueError("A polygon must have at least 3 points to calculate the area.")

        area = 0
        for i in range(n):
            x1, y1 = vertices[i].get_latitude(), vertices[i].get_longitude()
            x2, y2 = vertices[(i + 1) % n].get_latitude(), vertices[(i + 1) % n].get_longitude()
            area += x1 * y2 - x2 * y1

        return abs(area) / 2
    
    def size(self) -> int:
        """
        This method calculates the number of points of the Delimitation object.

        Returns: int

        Complexity: O(1)
        """
        return len(self.get_points())
    
    def add_point(self, point: Point):
        """
        This method adds a Point object to the Delimitation object.
        
        Arguments: point(Point)

        Complexity: O(1)
        """
        if type(point)!=Point:
            raise ValueError("Only Point objects can be added to the Delimitation")
        else: 
            if self.size()>0:
                self.__points.append(point)
                self.__third_last_point=self.__second_last_point
                self.__second_last_point=self.__last_point
                self.__last_point= point
            else:
                self.__points.append(point)
                self.__first_point=point
                self.__last_point= point
        
    def pop_point(self) -> Point:
        """
        This method removes thelast Point object that had been added to the
        Delimitation object.
        
        Returns: Point

        Complexity: O(1)
        """
        self.__second_last_point=self.__third_last_point
        self.__last_point=self.__second_last_point
        return self.__points.pop()

    def copy(self):
        """
        This method performs a deep copy of the Delimitation object.

        Returns: Delimitation

        Complexity: O(p), p being the number of points
        """
        d = Delimitation()
        points=self.get_points()
        for i in range(len(points)):
            d.add_point(points[i])
        return d
    
    def intersects(self, p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
        """
        This method determines whether the segment formed by the first two Point objects
        intersects the segment formed by the latter two Point objects.

        Arguments: p1(Point), p2(Point), p3(Point), p4(Point)

        Returns: bool

        Complexity: O(1)        
        """
        a1=(p2.get_longitude()-p1.get_longitude())
        a2=(p2.get_latitude()-p1.get_latitude())
        b1=(p4.get_longitude()-p3.get_longitude())
        b2=(p4.get_latitude()-p3.get_latitude())
        c1 = (a1 * p1.get_latitude() + b1 * p1.get_longitude())
        c2 = (a2 * p3.get_latitude() + b2 * p3.get_longitude())
        
        det= a1*b2-b1*a2

        if det==0:
            return False
        else:
            x = (b2 * c1 - b1 * c2) / det
            y = (a1 * c2 - a2 * c1) / det

            return (min(p1.get_longitude(), p2.get_longitude()) <= y <= max(p1.get_longitude(), p2.get_longitude()) and
            min(p1.get_latitude(), p2.get_latitude()) <= x <= max(p1.get_latitude(), p2.get_latitude()))  
    
    def crosses_delimitation(self, p1: Point, p2: Point) -> bool:
        """
        This method determines whether the segment formed by the two points intersects the
        delimitation. If the new segment ends in the first point of the delimitation,
        while not intersecting any other segment belonging to the delimitation, that this is not
        considered to be a crossing.
        
        Arguments: p1(Point), p2(Point)

        Returns: bool

        Complexity: O(p), p being the number of points
        """
        points=self.get_points()
        for i in range(1, len(points)):
            p3 = points[i - 1]
            p4 = points[i]

            if self.intersects(p3, p4, p1, p2):
                if p2 == self.get_first():
                    continue
                else:
                    return True
        return False
            
    def show(self, points: ValidPoint):
        """
        This method displays a window where all the stored points are displayed along
        with their identifier, and the complete delimitation formed by lines in between
        (possibly a subset of the) points is plotted.
        
        Arguments: points(ValidPoint)

        Complexity: O(p+v), p being the number of points in the Delimitation object
        and v being the number of points in the ValidPoint object.
        """
        coords = []
        for point in points.get_all_points():
            coords.append((point.get_latitude(), point.get_longitude(), point.get_id()))

        fig, ax = plt.subplots()

        if coords:
            latitudes, longitudes, point_ids = zip(*coords)
            ax.scatter(latitudes, longitudes, color='red')

            delim_latitudes = [point.get_latitude() for point in self.get_points()]
            delim_longitudes = [point.get_longitude() for point in self.get_points()]

            if len(delim_latitudes) > 1:
                ax.plot(delim_latitudes, delim_longitudes, color='blue', linewidth=1)

                first_point = self.get_first()
                last_point = self.get_last_two()[1]
                if first_point and last_point:
                    ax.plot(
                        [first_point.get_latitude(), last_point.get_latitude()],
                        [first_point.get_longitude(), last_point.get_longitude()],
                        color='blue', linewidth=1
                    )

            for i, point_id in enumerate(point_ids):
                ax.annotate(point_id, (latitudes[i], longitudes[i]), textcoords="offset points", xytext=(0, 5), ha='center')

        ax.set_title('Delimitation')
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Longitude')
        plt.show()

    def __eq__(self, del2) -> bool:
        """
        This method checks if the Delimitation objects that are being compared
        represent the same complete polygon.

        Arguments: del2(Delimitation)

        Returns: bool

        Complexity: O(p+m), p being the number of points in the first Delimitation 
        object and m being the number of points in the seconde Delimitation object.
        """
        points1=self.get_points()
        points2=del2.get_points()
        set1 = set(points1)
        set2 = set(points2)
        area1=self.get_area()
        area2=del2.get_area()
        return set1==set2 and area1==area2

    def __repr__(self) -> str:
        """
        This method returns a representation of the Delimitation object.

        Returns: str
        
        Complexity: O(p), p being the number of points in the Delimitation object
        """
        points = self.get_points()
        elements = ""
        for i in range(self.size()):
            if points[i]==self.get_last_two()[1]:
                elements += f"{points[i]}"
            else:
                elements += f"{points[i]} <- "
        return elements

