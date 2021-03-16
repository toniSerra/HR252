import copy
import math
import time
import matplotlib.pyplot as plt


def getpoints(n, s0, m):
    """
    given parameters s0, m return a list of n pseudorandom points
    :param n:   int
    :param s0:  int
    :param m:   int
    :return: list of points [px (int),py (int)]
    """
    newpoints = []
    while n > len(newpoints):
        s1 = (s0 ** 2) % m
        t1 = s1 % 2000 - 1000
        s2 = (s1 ** 2) % m
        t2 = s2 % 2000 - 1000
        newpoints.append((t1, t2))
        s0 = s2
    return newpoints


# 3 points [[ax, ay], [bx, by], [cx, cy]]
def areatriangle(a, b, c):
    # calcula l'area d'un triangles donada una llista amb 3 punts [[ax, ay], [bx, by], [cx, cy]]
    # Retorna el producte vectorial dels vectors AB i AC
    # producte vectorial = abs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay))/2
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2


# return True if P is in domain ABC
def outdomain(a, b, c, p):
    if p[0] <= min(a[0], b[0], c[0]) or p[0] >= max(a[0], b[0], c[0]) or p[1] <= min(a[1], b[1], c[1]) or p[1] >= max(a[1], b[1], c[1]):
        return True
    else:
        return False


# Boolean return True if P inside triangles ABC
def pointintriangle(a, b, c, p) -> bool:
    # comprovació si el punt p es troba dins el domini dels punts a, b i c
    # Si l'area del triangle ABC és igual als subtriangles PBC + APC + ABP, el punt P es troba a dins
    if areatriangle(a, b, c) == areatriangle(p, b, c) + areatriangle(a, p, c) + areatriangle(a, b, p):
        return True
    else:
        return False


# return angle's vector AB
def angle_vect(a, b):
    ang = math.atan2(b[1] - a[1], b[0] - a[0])
    if ang < 0:
        ang += 2 * math.pi
    return math.degrees(ang)


# forces angles to be in the rang(0-360)
def angle_norm(ang):
    ang %= 360
    if ang < 0:
        ang += 360
    return ang


# angles 3p
def anglestriangle(a, b, c):
    alfa_a = abs(angle_vect(a, c) - angle_vect(a, b))
    if alfa_a > 180:
        alfa_a = 360 - alfa_a
    alfa_b = abs(angle_vect(b, a) - angle_vect(b, c))
    if alfa_b > 180:
        alfa_b = 360 - alfa_b
    alfa_c = abs(angle_vect(c, b) - angle_vect(c, a))
    if alfa_c > 180:
        alfa_c = 360 - alfa_c
    if 180.0001 > alfa_a + alfa_b + alfa_c > 179.9999:
        return [alfa_a, alfa_b, alfa_c]
    else:
        print("ERROR! Algun triangle té uns angles que no sumen 180 +-10^-4", a,b, c)


def sorttrianglepoints(a, b, c):
    # given 3 indexpoints, sort them anticlockwise.
    center = [(pointlist[a][0] + pointlist[b][0] + pointlist[c][0])/3, (pointlist[a][1] + pointlist[b][1] + pointlist[c][1])/3]
    angle_lst = [angle_vect(center, pointlist[a]), angle_vect(center, pointlist[b]), angle_vect(center, pointlist[c])]
    angle_lst_sorted = sorted(angle_lst)
    point_lst = [a, b, c]
    point_lst_sorted = []

    for angle in angle_lst_sorted:
        point_lst_sorted.append(point_lst[angle_lst.index(angle)])

    return point_lst_sorted


# Return List of all valid triangles
def get_triangles(lst):
    # a, b and c are indexes of the list
    # check if any other given points are inside the triangle abc
    global areaMax
    global polygons_list
    global polygons_count
    global polyMax

    triangle_list = []
    triangles = 0

    el = len(lst)
    for a in range(0, el - 2):
        for b in range(a + 1, el - 1):
            for c in range(b + 1, el):
                triangle = True
                p = 0
                while triangle and p < el:
                    if p == a or p == b or p == c:
                        p += 1
                    elif outdomain(lst[a], lst[b], lst[c], lst[p]):
                        p += 1
                    else:
                        triangle = not pointintriangle(lst[a], lst[b], lst[c], lst[p])
                        p += 1

                if triangle:
                    triangle_indexes = sorttrianglepoints(a, b, c)
                    triangle_angles = anglestriangle(lst[triangle_indexes[0]], lst[triangle_indexes[1]],
                                                     lst[triangle_indexes[2]])
                    area = areatriangle(lst[a], lst[b], lst[c])
                    if area > areaMax:
                        areaMax = area
                        polyMax = [lst[a], lst[b], lst[c]]
                    triangles += 1
                    triangle_list.append([triangle_indexes, triangle_angles, area])

    polygons_count.append(triangles)
    polygons_count[0] += triangles  # 0 index refers to 'total'
    return triangle_list


def get_related_triangles(trilist):
    """
    Returns a paral·lel list to trianglelist with a list the indexes of the triangles that can be combined
    """
    reltriangles = []
    for t1 in range(len(trilist)):
        reltriangles.append([])

    for t1 in range(len(trilist)):
        tset = set(trilist[t1][0])
        for t2 in range(t1+1, len(trilist)):
            iset = set(trilist[t2][0])
            if len(tset.intersection(iset)) == 2:
                p = iset.difference(tset).pop()
                index = trilist[t2][0].index(p)
                ab = trilist[t2][0][index+1:]+trilist[t2][0][0:index] # common edge sorted anticlockwise relative to triangle i
                c = tset.difference(iset).pop()
                # pointP = pointlist[p]
                # pointA = pointlist[ab[0]]
                # pointB = pointlist[ab[1]]
                # pointC = pointlist[c]
                # print("points [A_B_C_P]  [", ab[0], ab[1], c, p,"]", pointA, pointB, pointC, pointP)
                angleAB = angle_norm(angle_vect(pointlist[ab[0]], pointlist[ab[1]]))
                # angleBC = angle_norm(angle_vect(pointB, pointC) - angleAB)
                # angleAC = angle_norm(angle_vect(pointA, pointC) - angleAB)
                # angleBP = angle_norm(angle_vect(pointB, pointP) - angleAB)
                # angleAP = angle_norm(angle_vect(pointA, pointP) - angleAB)
                if angle_norm(angle_vect(pointlist[ab[0]], pointlist[p])-angleAB) <\
                        angle_norm(angle_vect(pointlist[ab[0]], pointlist[c])-angleAB)-180\
                    and angle_norm(angle_vect(pointlist[ab[1]], pointlist[p])-angleAB) >\
                        angle_norm(angle_vect(pointlist[ab[1]], pointlist[c])-angleAB) - 180:
                    reltriangles[t1].append(t2)
                    reltriangles[t2].append(t1)
    return reltriangles


# Return True if 2 points match and evaluate the no matching point for checking if a new convex polygon can be created.
def triangleadjacent(pol, tri):
    pol_pt = list(copy.deepcopy(pol[0]))
    pol_an = list(copy.deepcopy(pol[1]))
    tri_pt = list(reversed(tri[0]))
    tri_an = list(reversed(tri[1]))
    for p in range(len(pol_pt)):
        for t in range(3):
            # if index's points matches:
            # [pol_pt[p], pol_pt[(p+1) % len(pol_pt)]] == [tri_pt[t], tri_pt[(t+1) % len(tri_pt)]]
            if pol_pt[p] == tri_pt[t] and pol_pt[(p+1) % len(pol_pt)] == tri_pt[(t+1) % 3]:
                angle1 = pol_an[p] + tri_an[t]
                angle2 = pol_an[(p + 1) % len(pol_an)] + tri_an[(t + 1) % len(tri_an)]
                # if both angles are acutes
                if angle1 < 180 and angle2 < 180:
                    pol_an[p] = angle1
                    pol_an[(p+1) % len(pol_an)] = angle2
                    pol_pt.insert(p + 1, tri_pt[t - 1])
                    pol_an.insert(p + 1, tri_an[t - 1])
                    # anticlockwise sorted points and angles begining with the minimum idex point.
                    ind = pol_pt.index(min(pol_pt))
                    pol_pt = pol_pt[ind:] + pol_pt[:ind]
                    pol_an = pol_an[ind:] + pol_an[:ind]
                    # print(pol_pt, pol_an)
                    return [pol_pt, pol_an, pol[2] + tri[2]]
    return False


def tmPrint(*args):
    if printTime:
        print(*args)


def dbPrint(*args):
    if printDebug:
        print(*args)


# CODE
printTime = True
printDebug = False

time0 = time.time()
polyMax = []            # biggest polygon
areaMax = 0             # biggest polygon's area
polygons_count = [0]    # cointains de number of polygons, total in index 0, 3 edges in idex 1 and so on.

# sample input parametres: 50, 290797, 50515093
pointlist = getpoints(30, 290797, 50515093)
triangles_containing_point = []  # point's paralel list with the triangle indexe's formed with each point

time1 = time.time()
tmPrint("t1_generació punts =", format(time1 - time0, '.2f'), "s")

triangles_list = get_triangles(pointlist)
time11 = time.time()
tmPrint("t11_generació triangles =", format(time11 - time1, '.2f'), "s")

triangles_related_triangles = get_related_triangles(triangles_list)
time12 = time.time()
tmPrint("t12_generació triangles relacionats =", format(time12 - time11, '.2f'), "s")

# copy for the first polygon generator iteration.
polygons_list = copy.deepcopy(triangles_list)
polygons_related_triangles = copy.deepcopy(triangles_related_triangles)

# index's triangles forming the polygon
# first iteration it would be just [0, 1, 2, 3,...] every triangle references his own position on triangle_list
# nexts iteration, every position will be filled with diferent combinations of triangles that could generate the polygon
polygons_inside_triangles = []
for p in range(len(polygons_list)):
    polygons_inside_triangles.append([p])

time2 = time.time()
tmPrint("t2_creació de llistes =", format(time2 - time12, '.2f'), "s")

edges = 4
while len(polygons_list) > 0:
    time31 = time.time()
    # The lists finishing with 'new' will replace the originals after every while's iteration
    polygons_list_new = []  # new poligons' list with n edges
    polygons_related_triangles_new = []  # polygons' related triangles' paral·lel list
    polygons_inside_triangles_new = []
    polygon_name_dict = {}  # new poligons name dict to check if the polygon is already created and its index
    polygon_count_new = 0

    for polygon_index in range(len(polygons_list)):
        polygon = polygons_list[polygon_index]
        related_triangles = polygons_related_triangles[polygon_index]  # triangles' index list
        dbPrint("checking polygon ", polygon_index, ":", polygon[0], "tri int", polygons_inside_triangles[polygon_index], "tri rel:", related_triangles)

        for related_triangle in range(len(related_triangles)):
            triangle_index = related_triangles[related_triangle]  # index of trianglelist
            triangle = triangles_list[triangle_index]  # triangle
            polygon_potencial = triangleadjacent(polygon, triangle)
            if polygon_potencial:
                dbPrint("\ttriangle {}: {} matches forming {}".format(triangle_index, triangle[0], polygon_potencial[0]))
                # join all de numbers to generate an unique String name
                ident = ""
                for i in polygon_potencial[0]:
                    ident += str(i) + " "

                if ident not in polygon_name_dict:
                    polygon_name_dict[ident] = polygon_count_new
                    polygon_count_new += 1
                    # Checking if it's the maximum polygon:
                    if polygon_potencial[-1] > areaMax:
                        areaMax = polygon_potencial[-1]
                        polyMax = polygon_potencial[0]
                    # Inside polygons:
                    # triangles_inside = polygons_inside_triangles[polygon_index]  # triangles adding to the list
                    polygons_inside_triangles_new.append(polygons_inside_triangles[polygon_index] + [triangle_index])
                    # Related polygons:
                    triangles_related = triangles_related_triangles[triangle_index]
                    triangles_rel_n_pol = set(polygons_related_triangles[polygon_index] +
                       triangles_related_triangles[triangle_index]).difference(set(polygons_inside_triangles_new[-1]))
                    dbPrint("\t\ttriangles related", polygons_related_triangles[polygon_index], "+",
                          triangles_related_triangles[triangle_index], "-",
                          polygons_inside_triangles_new[-1], "=",
                          list(triangles_rel_n_pol))
                    # Adding values to general lists
                    polygons_related_triangles_new.append(list(triangles_rel_n_pol))
                    polygons_list_new.append(polygon_potencial)
                    dbPrint("\t\tnew ident:", ident,
                          "\ttri ins:", polygons_inside_triangles_new[-1],
                          "\ttri rel", polygons_related_triangles_new[-1])
                else:
                    # despite being not unique polygon, inside triangles are added so as not have to check them again.
                    polygon_reference = polygon_name_dict[ident]            # index's polygon already existing
                    triangles_inside = polygons_inside_triangles[polygon_index] + [triangle_index]     # triangles adding to the list
                    polygons_inside_triangles_new[polygon_reference] = list(set(polygons_inside_triangles_new[polygon_reference]+triangles_inside))
                    dbPrint("\t\texistint ident:", ident,
                          "inner triangles to add:", triangles_inside)
                    # al segon quadrilater els poligons interiors haurien de ser 2,4,3,7
            else:
                x=1
                # afegir l'acció de borrar el triangle comprovat de la llista per no arrossegar-lo innecessaria.
                dbPrint("\ttriangle {}: {} doesn't match".format(triangle_index, triangle[0]))

    # Saving data for the next iteration.
    polygons_list = copy.deepcopy(polygons_list_new)
    polygons_related_triangles = list(polygons_related_triangles_new)
    polygons_inside_triangles = copy.deepcopy(polygons_inside_triangles_new)

    polygons_count.append(polygon_count_new)
    polygons_count[0] += polygon_count_new
    time32 = time.time()
    tmPrint("t3_generació ", edges, " =", format(time32 - time31, '.2f'), "s  ________________________________________")
    edges += 1


# stdout
print("%.1f" % areaMax)
print(polygons_count[0])

dbPrint("polygons:", polygons_count[1:], " total:", polygons_count[0])
dbPrint("Poligon maxim = ", polyMax)
tmPrint("temps d'execució =", format(time.time() - time0, '.2f'), "\bs")

# PLOTTING VARS
plot = 1                # 0,1
points = 1              # 0,1
pointlabels = 0         # 0,1
pointsize = 3           # 0-x
triangleall = 0         # 0,1

if plot:
    if points:
        index = 0
        for i in pointlist:
            plt.plot(i[0], i[1], marker='.', ms=pointsize, c='black')
            if pointlabels:
                plt.text(i[0] + 10, i[1] - 20, "P" + str(index))
            index += 1

    polyMaxX = []
    polyMaxY = []
    for p in range(len(polyMax)):
        polyMaxX.append(pointlist[polyMax[p]][0])
        polyMaxY.append(pointlist[polyMax[p]][1])
    polyMaxX.append(pointlist[polyMax[0]][0])
    polyMaxY.append(pointlist[polyMax[0]][1])
    plt.plot(polyMaxX,polyMaxY, c='red', linewidth=1, )

    plt.axis([-1000, 1000, -1000, 1000])
    plt.xlabel('eix X')
    plt.ylabel('eix Y')
    plt.show()
