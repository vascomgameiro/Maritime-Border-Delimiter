from typing import List

import numpy as np
from pythonds3.graphs import Graph

from daup2425p2 import Delimitation, Point, ValidPoints


class AproxOptimalDelimitation:
    """
    A class to calculate the approximate optimal delimitation of a set of points.
    """

    def __init__(self, valid_points: ValidPoints, max_distance: float):
        """
        Initializes the AproxOptimalDelimitation class with a set of valid points and a maximum distance.

        Args:
            valid_points (ValidPoints): A collection of valid points.
            max_distance (float): The maximum distance between connected points.

        Complexity: O(1)
        """
        self.__max_distance: float = max_distance
        self.__valid_points: ValidPoints = valid_points
        self.__graph: Graph = Graph()

    def build_graph(self):
        """
        Builds a graph connecting points within the maximum distance.

        Args:
            None

        Complexity: O(N^2)
        """
        points = self.__valid_points.get_all_points()

        for i, point1 in enumerate(points):
            for point2 in points[i + 1 :]:
                distance = point1.distance(point2)
                if distance <= self.__max_distance:
                    self.__graph.add_edge(point1, point2, weight=distance)

    def get_adjacency_list(self):
        """
        Generates an adjacency list where each point is connected to the points within the maximum distance.

        Points with no connections are skipped and not added to the adjacency list.

        Args:
            None

        Returns:
            dict: A dictionary where keys are Point objects, and values are lists of Point objects they are connected to.

        Complexity: O(N)
        """
        adjacency_list = {}
        invalid_points = set()

        for vertex in self.__graph:
            point: Point = vertex.get_key()  # Point object (this vertex)
            connections: List[Point] = [
                conn.get_key() for conn in vertex.get_neighbors() if conn.get_key() not in invalid_points
            ]  # Neighbour points

            if len(connections) == 1:
                invalid_points.add(point)
            elif len(connections) > 1:
                adjacency_list[point.get_id()] = connections

        return adjacency_list

    def find_delimitation(self) -> Delimitation:
        """
        Finds and returns the approximate optimal delimitation for the given set of points.

        The method builds a graph of valid points, prunes points with fewer than two valid connections,
        and then uses an angle-based traversal to approximate the optimal delimitation. The traversal
        continues until it loops back to the starting point.

        Complexity: O(N^2)

        Returns:
            Delimitation: The constructed delimitation connecting the points in an optimal manner.
        """

        self.build_graph()
        adjacency_data = self.get_adjacency_list()

        non_visited_points: List[Point] = [
            point for point in self.__valid_points.get_all_points() if point.get_id() in adjacency_data
        ]

        start_point: Point = min(non_visited_points, key=lambda p: p.get_longitude())
        delimitation = Delimitation()
        delimitation.add_point(start_point)

        current_point: Point = start_point
        reference_point = Point("ref", start_point.get_latitude() - 1, start_point.get_longitude() - 1)

        angle_dictionary = dict()
        for point in non_visited_points:
            angle = start_point.get_forward_angle(reference_point, point)
            angle_dictionary[point] = angle

        angle_dictionary.pop(start_point)
        while True:
            if delimitation.size() > 1:
                reference_point, current_point = delimitation.get_last_two()

            if not non_visited_points:
                return Delimitation()

            next_point = self._find_next_point_by_angle(
                current_point, reference_point, non_visited_points, adjacency_data
            )

            if next_point is None:
                non_visited_points.remove(current_point)
                delimitation.pop_point()

            else:
                delimitation.add_point(next_point)
                non_visited_points.remove(next_point)

                if next_point == start_point:
                    return delimitation

    def _find_next_point_by_angle(
        self, current_point: Point, reference_point: Point, non_visited_points: List[Point], adjacency_data: dict
    ):
        """
        Finds the next point by calculating the forward angle between the current point and each candidate.
        The point with the smallest angle is selected as the next point in the traversal.

        Args:
            current_point (Point): The current point in the traversal.
            reference_point (Point): The reference point used to calculate angles.
            non_visited_points (List[Point]): The list of candidate points to select from.

        Returns:
            Point: The next point in the traversal, or None if no valid point is found.

        Complexity:
            O(N)
        """
        min_angle = np.inf
        next_point = None

        candidate_points = [
            point
            for point in non_visited_points
            if (point in adjacency_data[current_point.get_id()] and point.get_id() in adjacency_data)
        ]

        if len(candidate_points) == 1:
            return candidate_points[0]
        if len(candidate_points) == 0:
            return None

        for point in candidate_points:
            angle = current_point.get_forward_angle(reference_point, point)
            if angle < min_angle:
                min_angle = angle
                next_point = point

        return next_point
