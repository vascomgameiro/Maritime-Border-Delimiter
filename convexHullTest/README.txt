fh_2 = FileHandler("test_.csv")
v_p2 = fh_2.read_points()
c = ConvexHull(v_p2)
c.find_delimitation()
