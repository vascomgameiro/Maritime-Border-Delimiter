# 🌊 Maritime Border Delimiter Project

## 📋 Overview
This project addresses the development of a Python tool to aid in the delimitation of maritime borders, particularly focusing on Portugal’s seabed territorial claims, in accordance with the United Nations Convention on the Law of the Sea (UNCLS). The goal is to create a polygonal delimitation that maximizes territory while following specific natural and political rules, simulating real-world negotiations with the UN.

The tool enables efficient handling of geological data points, drawing optimal or approximate delimitation boundaries, and producing maps for visualization.

## 🏆 Project Goals
1. **Efficient Delimitation**: Define a polygonal boundary that maximizes seabed territory while respecting distance constraints between points.
2. **Data Management**: Handle data from geological surveys and remote measurements, applying custom algorithms to define valid delimitation points.
3. **Visualization**: Automatically generate HTML maps to visualize seabed boundaries.

## 📐 Project Structure and Key Features
The project is divided into two parts, each with its specific focus:
- **Part 1**: Develop a command-line interface, read point data, and design Abstract Data Types (ADTs) for file handling, valid points, and point manipulation.
- **Part 2**: Implement algorithms for optimal and approximate delimitations, support different delimitation strategies, and extend ADTs to manage convex hulls and optimal paths.

## ⚙️ Commands and Usage
The tool operates through a command-line interface, processing a cycle of user-entered commands to perform tasks such as importing data, generating maps, and drawing delimitations.

### Main Commands
1. **`import_points <file_name>`**: Imports seabed point data from a CSV file and loads it into memory.
2. **`make_map`**: Generates an HTML map displaying all imported points and the latest calculated delimitation (if any), with a timestamped filename.
3. **`draw_hull`**: Calculates the convex hull that wraps all valid points, prints the boundary, its area, and computation time.
4. **`draw_optimal_delimitation <max_distance>`**: Finds a maximally large delimitation that respects a maximum distance between points.
5. **`draw_approx_delimitation <max_distance>`**: Efficiently approximates the optimal delimitation, focusing on time efficiency.
6. **`quit`**: Exits the program with a closing message.

## 📂 Abstract Data Types (ADTs)
The tool is designed with a set of specialized ADTs to handle key aspects of data manipulation and polygonal operations:
- **`FileHandler`**: Manages file interactions, reading point data, and creating HTML maps.
- **`ValidPoints`**: Stores and retrieves imported points with various retrieval and analysis functionalities.
- **`Point`**: Represents each seabed point, including geodesic distance calculations.
- **`Delimitation`**: Manages the points forming the delimitation, ensuring points are connected as a polygon.
- **Extended ADTs** (Part 2):
  - **`ConvexHull`**: Calculates the convex hull for a set of points.
  - **`OptimalDelimitation`** and **`ApproxDelimitation`**: Handle optimal and approximate delimitation calculations under specified distance constraints.

## 🧪 Analysis and Complexity
Each method and function within the ADTs includes a detailed docstring with complexity analysis, expressed in Big-O notation. This facilitates understanding of the tool's performance and ensures efficiency across large datasets.

