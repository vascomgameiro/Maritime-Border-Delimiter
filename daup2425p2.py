"""
This module contains classes and methods for processing geographical points and calculating distances.
"""

import math
import time
from datetime import datetime

import folium
import matplotlib.pyplot as plt
from geopy.distance import geodesic
from geopy.units import nautical
from pythonds3 import Graph, Vertex


class Point:
    """
    Represents a point on the globe with an identifier, latitude, and longitude.

    Internal Representation:,
    - __id (str): Unique identifier of the point,
    - __latitude (float): Latitude of the point,
    - __longitude (float): Longitude of the point,
    """

    def __init__(self, point_id: str, latitude: float, longitude: float) -> None:
        """
        Initializes a new Point object.

        Args:
            point_id (str): Identifier of the point.
            latitude (float): Latitude of the point.
            longitude (float): Longitude of the point.

        Complexity:
            O(1)
        """
        self.__id = point_id
        self.__latitude = latitude
        self.__longitude = longitude
        earth_radius = 6371009  # in meters
        lat_dist = math.pi * earth_radius / 180.0
        self.__y = nautical(meters=self.get_latitude() * lat_dist)  # type: ignore
        self.__x = nautical(meters=self.get_longitude() * lat_dist * math.cos(math.radians(self.get_latitude())))  # type: ignore

    def get_proj_x(self) -> float:
        """
        Returns the projected x-coordinate of the point.

        Returns:
            float: Projected x-coordinate of the point.

        Complexity:
            O(1)
        """
        return self.__x

    def get_proj_y(self) -> float:
        """
        Returns the projected y-coordinate of the point.

        Returns:
            float: Projected y-coordinate of the point.

        Complexity:
            O(1)
        """
        return self.__y

    def get_id(self) -> str:
        """
        Returns the identifier of the point.

        Returns:
            str: Identifier of the point.

        Complexity:
            O(1)
        """
        return self.__id

    def get_latitude(self) -> float:
        """
        Returns the latitude of the point.

        Returns:
            float: Latitude of the point.

        Complexity:
            O(1)
        """
        return self.__latitude

    def get_longitude(self) -> float:
        """
        Returns the longitude of the point.

        Returns:
            float: Longitude of the point.

        Complexity:
            O(1)
        """
        return self.__longitude

    def distance(self, other: "Point") -> float:
        """
        Calculates the geodesic distance between this point and another point.

        Args:
            other (Point): Another point to calculate the distance to.

        Returns:
            float: Distance in nautical miles between the two points.

        Complexity:
            O(1) - Uses geodesic calculation which is constant time.

        Note:
            The distance is returned in nautical miles.
        """
        coords_1 = (self.get_latitude(), self.get_longitude())
        coords_2 = (other.get_latitude(), other.get_longitude())
        distance = geodesic(coords_1, coords_2).nm

        return distance

    def in_vicinity_of(self, other: "Point", max_distance: int) -> bool:
        """
        Determines if this point is within a certain distance of another point.

        Args:
            other (Point): Another point to compare with.
            max_distance (int): Maximum distance in nautical miles.

        Returns:
            bool: True if within max_distance, False otherwise.

        Complexity:
            O(1)
        """
        return self.distance(other) <= max_distance

    def get_forward_angle(self, point_b: "Point", point_c: "Point") -> float:
        """
        Calculates the forward azimuth angle from this point to point_b,
        then to point_c.

        Args:
            point_b (Point): The second point in the sequence.
            point_c (Point): The third point in the sequence.

        Returns:
            float: The angle in degrees between the lines self->point_b and self->point_c.

        Complexity:
            O(1)

        Note:
            The angle is measured clockwise from the line point_b->self to self->point_c.
        """

        vector1 = (
            point_b.get_proj_x() - self.get_proj_x(),
            point_b.get_proj_y() - self.get_proj_y(),
        )

        vector2 = (
            point_c.get_proj_x() - self.get_proj_x(),
            point_c.get_proj_y() - self.get_proj_y(),
        )

        angle = math.degrees(math.atan2(vector1[1], vector1[0]) - math.atan2(vector2[1], vector2[0]))
        return (360 + angle) % 360

    def __repr__(self) -> str:
        """
        Returns the string representation of the Point object.

        Returns:
            str: String representation of the point.

        Complexity:
            O(1)
        """
        return f"Point {self.get_id()} ({self.get_latitude()}, {self.get_longitude()})"

    def __eq__(self, other: "Point") -> bool:
        """
        Checks if this point is equal to another point.

        Args:
            other (Point): The Point to compare with.

        Returns:
            bool: True if equal, False otherwise.

        Complexity:
            O(1)
        """
        if not isinstance(other, Point):
            return False
        return (
            self.get_id() == other.get_id()
            and self.get_latitude() == other.get_latitude()
            and self.get_longitude() == other.get_longitude()
        )

    def __hash__(self):
        """
        Returns the hash value of the Point object.

        The hash is based on the identifier, latitude and longitude of the point.

        Returns:
            int: The hash value of the point.

        Complexity:
            O(1)
        """
        return hash((self.__id, self.__latitude, self.__longitude))

    def __str__(self):
        """
        Returns the string representation of the Point object.

        Returns: str

        Complexity: O(1)
        """
        return f"({self.get_id()}, {self.get_latitude()}, {self.get_longitude()})"


class ValidPoints:
    """
    Represents a collection of valid points.

    Internal Representation:
    - __points (List[Point]): A list of Point objects.

    Methods:
        - get_size: Returns the number of points stored.
        - get_min_lat, get_max_lat, get_min_lon, get_max_lon: Return the min/max latitude/longitude.
        - get_all_points: Returns all the points stored.
        - get_points_vicinity: Returns the set of points within a given distance from a specified point.
    """

    def __init__(self, points: list):
        """
        Initializes the ValidPoints object with a list of Point objects.

        Args:
            points (List[Point]): A list of valid Point objects.

        Complexity:
            O(P) where P is the number of points in the list.
        """
        self.__points = points

    def get_size(self) -> int:
        """
        Returns the number of points stored.

        Returns:
            int: The number of points.

        Complexity:
            O(1)
        """
        return len(self.__points)

    def get_min_lat(self) -> float:
        """
        Returns the minimum latitude among the points stored.

        Returns:
            float: The minimum latitude.

        Complexity:
            O(P), where P is the number of points.
        """
        return min(point.get_latitude() for point in self.__points)

    def get_max_lat(self) -> float:
        """
        Returns the maximum latitude among the points stored.

        Returns:
            float: The maximum latitude.

        Complexity:
            O(P), where P is the number of points.
        """
        return max(point.get_latitude() for point in self.__points)

    def get_min_lon(self) -> float:
        """
        Returns the minimum longitude among the points stored.

        Returns:
            float: The minimum longitude.

        Complexity:
            O(P), where P is the number of points.
        """
        return min(point.get_longitude() for point in self.__points)

    def get_max_lon(self) -> float:
        """
        Returns the maximum longitude among the points stored.

        Returns:
            float: The maximum longitude.

        Complexity:
            O(P), where P is the number of points.
        """
        return max(point.get_longitude() for point in self.__points)

    def get_all_points(self) -> list:
        """
        Returns the list of all points stored, sorted by their identifier.

        Returns:
            List[Point]: A list of all Point objects.

        Complexity:
            O(P log P) where P is the number of points (due to sorting).
        """
        return sorted(self.__points, key=lambda point: point.get_id())

    def get_points_vicinity(self, ref_point: "Point", max_distance: int) -> set:
        """
        Returns a set of points that are within a given maximum distance from the reference point.

        Args:
            ref_point (Point): The reference point.
            max_distance (int): The maximum distance in nautical miles.

        Returns:
            Set[Point]: A set of points that are within the max_distance from the reference point.

        Complexity:
            O(P), where P is the number of points (since we must check each point).
        """
        return {point for point in self.__points if point.in_vicinity_of(ref_point, max_distance)}

    def __repr__(self) -> str:
        """
        Returns the string representation of the ValidPoints object, showing all points.

        Returns:
            str: A string representation of the points stored.

        Complexity:
            O(P), where P is the number of points.
        """
        return f"ValidPoints({self.__points})"


class Delimitation:
    """
    The Delimitation ADT models a polygon that can be constructed by sequentially adding
    points.

    Attributes:
        __points (list): A List that stores all the points in the order they are added.

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
        show(ValidPoints): Displays the current Delimitation and points using a plotting tool.
    """

    def __init__(self):
        """
        This method initializes an empty Delimitation object.

        Complexity: O(1)
        """
        self.__points = []

    def get_points(self) -> list:
        """
        This method retrieves the Points that belong to the Delimitation object by the order
        they were added.

        Returns: list

        Complexity: O(1)
        """
        return self.__points

    def get_first(self) -> Point | None:
        """
        This method retrieves the first point that was added to the Delimitation object.

        Returns: Point

        Complexity: O(1)
        """
        return self.__points[0] if self.size() > 0 else None

    def get_last_two(self) -> tuple:
        """
        This method retrieves the last two Point objects that have been added to the
        Delimitation object,the rightmost being the last one added.

        Returns: tuple(Point, Point)

        Complexity: O(1)
        """
        if self.size() == 0:
            return (None, None)
        elif self.size() == 1:
            return (None, self.__points[-1])
        return (self.__points[-2], self.__points[-1])

    def get_area(self) -> float:
        """
        This method calculates the area of the polygon that is represented. If
        the Delimitation object has less than 3 points or collinear points,
        the method raises a ValueError.

        Returns: float

        Complexity: O(p), p being the number of points
        """
        vertices = self.get_points()
        n = len(vertices)

        if n < 3:
            raise ValueError("A polygon must have at least 3 points to calculate the area.")

        if self.__collinearity_test(vertices):
            raise ValueError("It is impossible to form a polygon with collinear points.")

        area = 0
        for i in range(n):
            x1, y1 = vertices[i].get_proj_x(), vertices[i].get_proj_y()
            x2, y2 = vertices[(i + 1) % n].get_proj_x(), vertices[(i + 1) % n].get_proj_y()
            area += x1 * y2 - x2 * y1

        return abs(area) / 2

    @staticmethod
    def __collinearity_test(points: list) -> bool:
        """
        This static method checks if a list of Point objects is collinear.

        Arguments: points(list)

        Returns: bool

        Complexity: O(p), p being the number of points
        """
        if len(points) < 3:
            return True

        x1, y1 = points[0].get_latitude(), points[0].get_longitude()
        x2, y2 = points[1].get_latitude(), points[1].get_longitude()

        for i in range(2, len(points)):
            x3, y3 = points[i].get_latitude(), points[i].get_longitude()
            if (x2 - x1) * (y3 - y1) != (y2 - y1) * (x3 - x1):
                return False

        return True

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
        if not isinstance(point, Point):
            raise ValueError("Only Point objects can be added to the Delimitation object")
        if point not in self.get_points():
            self.__points.append(point)

    def pop_point(self) -> Point:
        """
        This method removes the last Point object that had been added to the
        Delimitation object.

        Returns: Point

        Complexity: O(1)
        """
        if self.size() == 0:
            raise IndexError("Cannot remove point from an empty delimitation.")
        return self.__points.pop()

    def copy(self) -> "Delimitation":
        """
        This method performs a deep copy of the Delimitation object.

        Returns: Delimitation

        Complexity: O(p), p being the number of points
        """
        delimitation_copy = Delimitation()
        points = self.get_points()

        for i in range(self.size()):
            point_copy = Point(points[i].get_id(), points[i].get_latitude(), points[i].get_longitude())
            delimitation_copy.add_point(point_copy)

        return delimitation_copy

    @staticmethod
    def __min_max_lat_lon(p1: Point, p2: Point) -> tuple:
        """
        Helper function to return min and max latitudes and longitudes for a pair of points.

        This function takes two Point objects as arguments and returns a tuple of two tuples.
        The first tuple contains the minimum and maximum latitude, and the second tuple
        contains the minimum and maximum longitude.

        Args:
            p1 (Point): The first point.
            p2 (Point): The second point.

        Returns:
            tuple: A tuple of two tuples containing the min and max latitude, and min and max longitude.
        """
        return (
            (min(p1.get_latitude(), p2.get_latitude()), max(p1.get_latitude(), p2.get_latitude())),
            (
                min(p1.get_longitude(), p2.get_longitude()),
                max(p1.get_longitude(), p2.get_longitude()),
            ),
        )

    def __is_within_segment(self, p1: Point, p2: Point, intersect_x: float, intersect_y: float) -> bool:
        """
        This static method determines whether a point is within a segment formed by two other points.

        Arguments: p1(Point), p2(Point), p3(Point)

        Returns: bool

        Complexity: O(1)
        """
        (x1_min, x1_max), (y1_min, y1_max) = self.__min_max_lat_lon(p1, p2)
        return (x1_min <= intersect_x <= x1_max) and (y1_min <= intersect_y <= y1_max)

    def intersects(self, p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
        """
        This method determines whether the segment formed by the first two Point objects
        intersects the segment formed by the latter two Point objects.

        Arguments: p1(Point), p2(Point), p3(Point), p4(Point)

        Returns: bool

        Complexity: O(1)
        """
        if not all(isinstance(p, Point) for p in [p1, p2, p3, p4]):
            raise ValueError("The four points must be of type 'Point'.")

        intersection = self.__intersection_point(p1, p2, p3, p4)

        if intersection is None:
            return False

        intersect_x, intersect_y = intersection.get_latitude(), intersection.get_longitude()

        is_within_first_segment = self.__is_within_segment(p1, p2, intersect_x, intersect_y)
        is_within_second_segment = self.__is_within_segment(p3, p4, intersect_x, intersect_y)

        return is_within_first_segment and is_within_second_segment

    @staticmethod
    def __intersection_point(p1: Point, p2: Point, p3: Point, p4: Point) -> Point | None:
        """
        This method determines the intersection point, if it exists, between the segment formed by the
        first two Point objects and the segment formed by the latter two Point objects.

        Arguments: p1(Point), p2(Point), p3(Point), p4(Point)

        Returns: Point

        Complexity: O(1)
        """
        a1 = p2.get_longitude() - p1.get_longitude()
        b1 = p1.get_latitude() - p2.get_latitude()
        c1 = a1 * (p1.get_latitude()) + b1 * (p1.get_longitude())

        a2 = p4.get_longitude() - p3.get_longitude()
        b2 = p3.get_latitude() - p4.get_latitude()
        c2 = a2 * (p3.get_latitude()) + b2 * (p3.get_longitude())

        det = a1 * b2 - a2 * b1
        if abs(det) < 1e-9:
            return None

        x = (b2 * c1 - b1 * c2) / det
        y = (a1 * c2 - a2 * c1) / det
        return Point("intersection", x, y)

    def crosses_delimitation(self, p1: Point, p2: Point) -> bool:
        """
        This method determines whether the segment formed by the two points intersects the
        delimitation. If the new segment ends in the first point of the delimitation,
        while not intersecting any other segment belonging to the delimitation, this is not
        considered to be a crossing.

        Arguments: p1(Point), p2(Point)

        Returns: bool

        Complexity: O(p), p being the number of points
        """
        if not all(isinstance(p, Point) for p in [p1, p2]):
            raise ValueError("Both points must be of type 'Point'.")

        if self.size() == 0:
            raise ValueError("Delimitation object is empty.")

        points = self.get_points()
        for i in range(1, len(points)):
            p3 = points[i - 1]
            p4 = points[i]

            if self.intersects(p1=p3, p2=p4, p3=p1, p4=p2):
                intersection = self.__intersection_point(p1=p3, p2=p4, p3=p1, p4=p2)
                intersect_x, intersect_y = intersection.get_latitude(), intersection.get_longitude()  # type: ignore
                first_point = self.get_first()
                last_point = self.get_last_two()[1]

                if (intersect_x == first_point.get_latitude() and intersect_y == first_point.get_longitude()) or (  # type: ignore
                    intersect_x == last_point.get_latitude() and intersect_y == last_point.get_longitude()  # type: ignore
                ):
                    continue
                return True

        return False

    def show(self, points: ValidPoints):
        """
        This method displays a window where all the stored points are displayed along
        with their identifier, and the complete delimitation formed by lines in between
        (possibly a subset of the) points is plotted.

        Arguments: points(ValidPoints)

        Complexity: O(p+v), p being the number of points in the Delimitation object
        and v being the number of points in the ValidPoints object.
        """
        if not isinstance(points, ValidPoints):
            raise ValueError("The stored points must be of type 'ValidPoints'.")

        coords = []
        for point in points.get_all_points():
            coords.append((point.get_proj_x(), point.get_proj_y(), point.get_id()))

        _, ax = plt.subplots()

        if coords:
            latitudes, longitudes, point_ids = zip(*coords)
            ax.scatter(latitudes, longitudes, color="red")

            for i, point_id in enumerate(point_ids):
                ax.annotate(
                    point_id,
                    (latitudes[i], longitudes[i]),
                    textcoords="offset points",
                    xytext=(0, 5),
                    ha="center",
                )

        if self.size() > 1:
            delim_latitudes = [point.get_proj_x() for point in self.get_points()]
            delim_longitudes = [point.get_proj_y() for point in self.get_points()]
            ax.plot(delim_latitudes, delim_longitudes, color="blue", linewidth=1)

            if len(delim_latitudes) > 2:
                first_point = self.get_first()
                last_point = self.get_last_two()[1]
                ax.plot(
                    [first_point.get_proj_x(), last_point.get_proj_x()],  # type: ignore
                    [first_point.get_proj_y(), last_point.get_proj_y()],  # type: ignore
                    color="blue",
                    linewidth=1,
                )

        ax.set_title("Delimitation")
        ax.set_xlabel("Latitude")
        ax.set_ylabel("Longitude")
        plt.show()

    def __eq__(self, other: "Delimitation") -> bool:
        """
        This method checks if the Delimitation objects that are being compared
        represent the same complete polygon.

        Arguments: other(Delimitation)

        Returns: bool

        Complexity: O(p^2), p being the number of points in the Delimitation object
        """
        if not isinstance(other, Delimitation):
            return False

        points1 = self.get_points()
        points2 = other.get_points()
        set1 = set(points1)
        set2 = set(points2)
        if set1 == set2:
            return self.__check_rotations(points1, points2)

        return False

    @staticmethod
    def __check_rotations(list1: list, list2: list) -> bool:
        """
        This static method checks if list2 is a rotation of list1, considering both
        clockwise and counterclockwise rotations.

        Arguments: list1 (list), list2 (list)

        Returns: bool

        Complexity: O(p^2), p being the number of points in the lists
        """
        double_list1 = list1 + list1
        for i in range(len(list1)):
            if double_list1[i : i + len(list1)] == list2:
                return True

        reverse_list2 = list2[::-1]
        for i in range(len(list1)):
            if double_list1[i : i + len(list1)] == reverse_list2:
                return True

        return False

    def is_valid_delimitation(self, valid_points: ValidPoints, distance: int) -> bool:
        """
        Validates if the current delimitation formed by the points of the object is 
        valid based on a set of conditions.

        Parameters: valid_points(ValidPoints), distance(int)

        Returns: bool

        Complexity: O(p^2), p being being the number of points in the Delimitation object
        """
        del_points = self.get_points()
        v_p_points = valid_points.get_all_points()
        if set(del_points).issubset(set(v_p_points)):
            if self.get_first().distance(self.get_last_two()[1]) > distance:
                return False
            for i in range(self.size() - 1):
                if del_points[i].distance(del_points[i+1]) <= distance:
                    continue
                else:
                    return False
            for i in range(len(del_points)):
                p1 = del_points[i]
                p2 = del_points[(i + 1) % len(del_points)] 

                for j in range(len(del_points)):
                    if j == i or (j + 1) % len(del_points) == i or (j == (i + 1) % len(del_points)):
                        continue 
                    
                    p3 = del_points[j]
                    p4 = del_points[(j + 1) % len(del_points)]  
                    
                    if self.intersects(p1, p2, p3, p4):
                        return False  

            return True
        return False

    def __repr__(self) -> str:
        """
        This method returns a representation of the Delimitation object.

        Return: str

        Complexity: O(p), p being the number of points in the Delimitation object
        """
        points_repr = ", ".join([repr(point) for point in self.get_points()])
        return f"Delimitation([{points_repr}])"

    def __str__(self):
        """
        Returns a string representation of the object

        Return: str

        Complexity: O(p), p being the number of points in the Delimitation object
        """
        points_repr = ", ".join([repr(point) for point in self.get_points()])
        return f"[{points_repr}]"


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
        val = (p2.get_longitude() - p1.get_longitude()) * (p3.get_latitude() - p2.get_latitude()) - (
            p2.get_latitude() - p1.get_latitude()
        ) * (p3.get_longitude() - p2.get_longitude())

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
        sorted_points = sorted(dictionary1, key=dictionary1.get)  # type: ignore

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

    def __dfs(self, current_vertex: Vertex, delimitation: Delimitation, area_registry: dict, visited: set):
        """
        Performs a Depth-First Search (DFS) traversal starting from the given vertex, building
        possible delimitations and calculating their areas. The valid delimitations are stored
        in `area_registry`.

        Arguments: current_vertex: The vertex from which DFS is initiated (a Vertex object).
                   delimitation: A Delimitation object that tracks the sequence of points in the current path.
                   area_registry: A Dictionary to store valid delimitation areas mapped to the corresponding Delimitation object.
                   visited: A Set to track visited vertices during the DFS traversal.

        Complexity: O(v+e+v*p), v being the number of vertices, e being the number of edges and p being the number of
        points in the delimitation object
        """
        visited.add(current_vertex)
        if current_vertex.get_neighbors() is not None:
            for next_vertex in current_vertex.get_neighbors():
                if next_vertex not in visited:
                    delimitation.add_point(next_vertex.get_key())
                    if delimitation.size() >= 3:
                        try:
                            total_area = delimitation.get_area()
                            area_registry[total_area] = delimitation.copy()
                        except ValueError:
                            continue
                    self.__dfs(next_vertex, delimitation, area_registry, visited)
                    if delimitation.size() > 0:
                        delimitation.pop_point()
            visited.remove(current_vertex)

    def find_delimitation(self) -> Delimitation:
        """
        Finds and returns the delimitation with the maximum area by building a graph of points and performing DFS.

        The method constructs a graph where points are vertices and edges represent the valid connections
        between points based on the distance threshold. DFS is performed on each point, and the valid
        delimitations are stored. The delimitation with the maximum area is returned.

        Returns: Delimitation object

        Complexity: O(p^2(1+vd)), p being the number of vertices in the graph and vd being the number of valid
        Delimitation objects

        """
        area_registry = dict()
        points = self.__valid_points.get_all_points()
        g = Graph()
        for i in range(self.__valid_points.get_size()):
            neigbors = self.__valid_points.get_points_vicinity(points[i], self.__distance)
            for neighbor in neigbors:
                if points[i] != neighbor:
                    g.add_edge(points[i], neighbor, points[i].distance(neighbor))

        for vertex_key in g.get_vertices():
            visited = set()
            d = Delimitation()
            d.add_point(vertex_key)  # type: ignore
            self.__dfs(g.get_vertex(vertex_key), d, area_registry, visited)  # type: ignore

        valid_area_registry = {
            area: delimitation
            for area, delimitation in area_registry.items()
            if delimitation.is_valid_delimitation(self.__valid_points, self.__distance)
        }
        maximum = max(valid_area_registry)
        return area_registry[maximum]


class FileHandler:
    """
    The FileHandler ADT is responsible for managing file interactions within the program.
    It handles reading a specified file containing data about points and generating HTML maps
    using those points.

    Attributes:
        __filename (str): The name of the file from which point data will be read.
        __points (ValidPoints): The ValidPoints object made of the file points

    Methods:
        __init__(str): Initializes a new FileHandler object with the specified file name.
        The file handler is automatically returned at instantiation.

        read_points(): Reads data about points from the specified file and
        creates a ValidPoints object to store the Point objects. Returns the ValidPoints object.

        make_map(Delimitation): Creates an HTML file containing a map
        that represents the points stored in the ValidPoints object as circles,
        and the provided Delimitation object as a line. Uses the folium module to generate the map.

        __repr__(): Returns a string representation of the FileHandler object, including
        the file name and the points that have been imported.
    """

    def __init__(self, filename: str):
        """
        This method initializes a FileHandler object.

        Arguments: filename(str)

        Complexity: O(1)
        """
        self.__filename = filename
        self.__points = self.read_points()

    def read_points(self) -> ValidPoints:
        """
        This method reads point data from a file and stores it in a ValidPoints object.

        Returns: ValidPoints

        Complexity: O(p), p being the number of points stored in the file
        """
        valid_points = []
        open_file = open(self.__filename, "r", encoding="utf-8")
        if open_file:
            for line in open_file:
                items = line.split()
                valid_points.append(Point(items[0], float(items[1]), float(items[2])))
        else:
            print(f"Error opening file {self.__filename}")

        open_file.close()
        return ValidPoints(valid_points)

    def make_map(self, delimitation=None):
        """
        This method creates a HTML file containing a map that represents the points stored
        in the file handler's ValidPoints object as circles, and the delimitation as
        a line.

        Arguments: delimitation(Delimitation)

        Complexity: O(plog(p)+d), p being the number of points of the ValidPoints object and d
        being the number of points of the Delimitation object
        """
        folium_map = folium.Map()
        points = self.__points.get_all_points()
        coords = [(point.get_latitude(), point.get_longitude()) for point in points]

        min_lat = self.__points.get_min_lat()
        min_lon = self.__points.get_min_lon()

        max_lat = self.__points.get_max_lat()
        max_lon = self.__points.get_max_lon()

        folium_map.fit_bounds([[min_lat, min_lon], [max_lat, max_lon]])

        for latitude, longitude in coords:
            folium.CircleMarker(location=[latitude, longitude], radius=2, weight=5).add_to(folium_map)

        if delimitation:
            if delimitation.size() > 0:
                del_points = delimitation.get_points()
                del_coordinates = [(del_point.get_latitude(), del_point.get_longitude()) for del_point in del_points]
                del_coordinates.append(
                    (
                        delimitation.get_first().get_latitude(),
                        delimitation.get_first().get_longitude(),
                    )
                )
                folium.PolyLine(locations=del_coordinates, color="red", weight=2.5).add_to(folium_map)

        self.__filename = f"map_{time.time()}.html"
        folium_map.save(self.__filename)

    def __repr__(self) -> str:
        """
        This method returns a representation of the FileHandler object.

        Returns: str

        Complexity: O(p), p being the number of points in the FileHandler object
        """
        return f"Filename: {self.__filename}; {self.__points}"


def main():
    """
    Main function of the program, which will continuously prompt the user for input commands.
    Commands available are:
    - import_points <file_name>
    - draw_hull
    - draw_optimal_delimitation <max_distance>
    - draw_approx_delimitation <max_distance>
    - make_map
    - quit
    """
    file_handler = None
    delimitation = None

    while True:
        command = input("?> ").strip().split()
        if len(command) == 0:
            continue
        cmd = command[0]

        if cmd == "import_points":
            if len(command) < 2:
                print("Error: missing file name.")
                continue
            file_name = command[1]

            try:
                start_time = time.time()
                file_handler = FileHandler(file_name)
                num_points = file_handler.read_points().get_size()
                elapsed_time = (time.time() - start_time) * 1000
                print(f"Size: {num_points}")
                print(f"Time: {elapsed_time:.2f}")
            except Exception as e:
                print(f"Error: {e}")
                continue

        elif cmd == "draw_hull":
            if file_handler is None:
                print("Error: No points imported. Please use 'import_points' first.")
                continue
            try:
                start_time = time.time()
                delimitation = ConvexHull(file_handler.read_points()).find_delimitation()
                area = delimitation.get_area()
                elapsed_time = (time.time() - start_time) * 1000

                print("Type: Convex Hull")
                print(f"Line: {delimitation}")
                print(f"Area: {area:.2f}")
                print(f"Time: {elapsed_time:.2f}")
            except Exception as e:
                print(f"Error: {e}")

        elif cmd == "draw_optimal_delimitation":
            if len(command) < 2:
                print("Error: missing max distance argument.")
                continue
            if file_handler is None:
                print("Error: No points imported. Please use 'import_points' first.")
                continue
            try:
                max_distance = int(command[1])
                start_time = time.time()
                delimitation = OptimalDelimitation(file_handler.read_points(), max_distance).find_delimitation()
                elapsed_time = (time.time() - start_time) * 1000
                area = delimitation.get_area()
                print("Type: Optimal Delimitation")
                print(f"Line: {delimitation}")
                print(f"Area: {area:.2f}")
                print(f"Time: {elapsed_time:.2f}")
            except Exception as e:
                print(f"Error: {e}")
        # TODO: Implement approximate delimitation
        # elif cmd == "draw_approx_delimitation":
        #     if len(command) < 2:
        #         print("Error: missing max distance argument.")
        #         continue
        #     if file_handler is None:
        #         print("Error: No points imported. Please use 'import_points' first.")
        #         continue
        #     try:
        #         max_distance = int(command[1])
        #         start_time = time.time()
        #         delimitation = ApproximateDelimitation(file_handler.read_points(), max_distance)
        #         elapsed_time = (time.time() - start_time) * 1000
        #         area = delimitation.calculate_area()
        #         print("Type: Approximate Delimitation")
        #         print(f"Line: {delimitation}")
        #         print(f"Area: {area:.2f}")
        #         print(f"Time: {elapsed_time:.2f}")
        #     except Exception as e:
        #         print(f"Error: {e}")

        elif cmd == "make_map":
            if file_handler is None or delimitation is None:
                print("Error: No points or delimitation available.")
                continue
            try:
                start_time = time.time()
                timestamp = datetime.now().timestamp()
                file_name = f"map_{timestamp}.html"
                file_handler.make_map(delimitation)
                elapsed_time = (time.time() - start_time) * 1000
                print(f"Map: {file_name}")
                print(f"Time: {elapsed_time:.2f}")
            except Exception as e:
                print(f"Error: {e}")

        elif cmd == "quit":
            print("Done!")
            break

        else:
            print(
                f"Unknown command: {cmd}. Valid Commands: import_points, draw_hull, draw_optimal_delimitation <max_distance>, draw_approx_delimitation <max_distance>, make_map, quit"
            )


if __name__ == "__main__":
    main()