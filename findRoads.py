import overpy, os, requests, time, random, json, statistics
from datetime import datetime
import streamlit as st
import pandas as pd
from math import radians, degrees, sin, cos, sqrt, atan2
from io import StringIO

debug = False
metric = True

#test/debug dataset controls
#origin
home = False
willow = False
outlets = False
interchange = False

outFolder = "debugData/"
fName = outFolder + "debugRoads" #base file name with no extension
jsonWrite = True #write method
if debug and not os.path.exists(outFolder):
    os.makedirs(outFolder)

class ThreeDimSegment:
    start = int
    end = int
    elevation = []

    def __init__(self, start, end, elevCoords):
        self.start = start
        self.end = end
        self.elevation = elevCoords

def calculate_bounding_box(latitude, longitude, radius):
    """
    Calculate the bounding box given a pair of coordinates (latitude and longitude)
    and a radius in miles.

    Args:
        latitude (float): Latitude of the center point (in degrees).
        longitude (float): Longitude of the center point (in degrees).
        radius (float): Radius in miles.

    Returns:
        tuple: A tuple containing (lowLatitude, leftLongitude, highLatitude, rightLongitude)
               representing the bounding box.
    """
    # Earth's radius in miles
    if metric:
        earthRadius = 6371.0 # kilometers
    else:
        earthRadius = 3958.8  # miles

    # Convert latitude and longitude from degrees to radians
    latRad = radians(latitude)
    lonRad = radians(longitude)

    # Convert radius from miles to radians
    radiusRad = radius / earthRadius

    # Calculate the minimum and maximum latitude
    lowLat = degrees(latRad - radiusRad)
    highLat = degrees(latRad + radiusRad)

    # Calculate the minimum and maximum longitude
    deltaLon = sqrt(radiusRad**2 - (sin(latRad) * sin(radiusRad))**2)
    leftLon = degrees(lonRad - deltaLon)
    rightLon = degrees(lonRad + deltaLon)

    return lowLat, leftLon, highLat, rightLon

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Convert decimal degrees to radians
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    # Radius of the Earth in feet
    if metric:
        R = 6371000.0
    else:
        R = 20902231.84  # feet
    distance = R * c
    return distance

def combineSegments(segmentList):
    combinedList = []
    combinedList.append(segmentList[0])
    del (segmentList[0])

    while (len(segmentList) > 0):
        minDistanceToEndpoint = float('inf')
        minToEndpointSegment = -1
        minDistanceToStartpoint = float('inf')
        minToStartpointSegment = -1
        toBeAttachedTo = -1
        minDistance = float('inf')
        rightLeft = True  #Right is true, left is false
        reverseEndpoint = False
        reverseStartpoint = False
        for j in range(0, len(combinedList)):
            #compare segmentList[0] startpoint to all endpoints of segmentList segments
            #save closest segment endpoint to startpoint && distance to startpoint
            #compare segmentList[0] endpoint to all startpoints of segmentList segments
            #save closest segment startpoint to endpoint && distance to endpoint
            #attach segment with minimum distance

            for i in range(0, len(segmentList)):
                distance = haversine(segmentList[i][0][0], segmentList[i][0][1], combinedList[j][-1][0], combinedList[j][-1][1])
                if distance < minDistanceToEndpoint:
                    #This uncombined segment's startpoint is closest to the current combined segment's endpoint
                    minDistanceToEndpoint = distance
                    minToEndpointSegment = i
                    reverseEndpoint = False

                distance = haversine(segmentList[i][-1][0], segmentList[i][-1][1], combinedList[j][-1][0], combinedList[j][-1][1])
                if distance < minDistanceToEndpoint:
                    #This uncombined segment's endpoint is closest to the current combined segment's endpoint
                    minDistanceToEndpoint = distance
                    minToEndpointSegment = i
                    reverseEndpoint = True

                distance = haversine(segmentList[i][-1][0], segmentList[i][-1][1], combinedList[j][0][0], combinedList[j][0][1])
                if distance < minDistanceToStartpoint:
                    #This uncombined segment's endpoint is closest to the current combined segment's startpoint
                    minDistanceToStartpoint = distance
                    minToStartpointSegment = i
                    reverseStartpoint = False

                distance = haversine(segmentList[i][0][0], segmentList[i][0][1], combinedList[j][0][0], combinedList[j][0][1])
                if distance < minDistanceToStartpoint:
                    #This uncombined segment's startpoint is closest to the current combined segment's startpoint
                    minDistanceToStartpoint = distance
                    minToStartpointSegment = i
                    reverseStartpoint = True

            if minDistanceToEndpoint < minDistanceToStartpoint:
                minDistance = minDistanceToEndpoint
                toBeAttachedTo = j
                rightLeft = True

            if minDistanceToStartpoint < minDistance:
                minDistance = minDistanceToStartpoint
                toBeAttachedTo = j
                rightLeft = False


        if rightLeft:
            if reverseEndpoint:
                segmentList[minToEndpointSegment].reverse()
            combinedList.insert(toBeAttachedTo + 1, segmentList[minToEndpointSegment])
            del (segmentList[minToEndpointSegment])

        elif not rightLeft:
            if reverseStartpoint:
                segmentList[minToStartpointSegment].reverse()
            if toBeAttachedTo > 0:
                combinedList.insert(toBeAttachedTo - 1,
                                    segmentList[minToStartpointSegment])
            else:
                combinedList.insert(toBeAttachedTo,
                                    segmentList[minToStartpointSegment])
            del (segmentList[minToStartpointSegment])

    finalList = []
    for seg in combinedList:
        for pair in seg:
            finalList.append(pair)

    return finalList

def getRoads(lowLat, leftLong, highLat, rightLong, pBar, p):
    api = overpy.Overpass()

    with st.spinner("Getting coordinates of all roads in bounding box"):
        tertiaryResult = api.query("""
            way(%s,%s,%s,%s) ["highway" = "tertiary"];
            (._;>;);
            out body;
            """ % (str(lowLat), str(leftLong), str(highLat), str(rightLong),))
        p += 1
        pBar.progress(p, "Getting coordinates...")
        motorwayResult = api.query("""
            way(%s,%s,%s,%s) ["highway" = "motorway"];
            (._;>;);
            out body;
            """ % (str(lowLat), str(leftLong), str(highLat), str(rightLong),))
        p += 1
        pBar.progress(p, "Getting coordinates...")
        secondaryResult = api.query("""
            way(%s,%s,%s,%s) ["highway" = "secondary"];
            (._;>;);
            out body;
            """ % (str(lowLat), str(leftLong), str(highLat), str(rightLong),))
        p += 1
        pBar.progress(p, "Getting coordinates...")
        primaryResult = api.query("""
            way(%s,%s,%s,%s) ["highway" = "primary"];
            (._;>;);
            out body;
            """ % (str(lowLat), str(leftLong), str(highLat), str(rightLong),))
        p += 1
        pBar.progress(p, "Sorting coordinates...")
    
    st.success("Retrieved road coordinates", icon='✔')
        

    
    results = []
    results.append(tertiaryResult)
    results.append(motorwayResult)
    results.append(secondaryResult)
    results.append(primaryResult)

    if debug:
        with open(fName + "AllSegments.txt", "w") as f:
            for result in results:
                for way in result.ways:
                    f.write("Name: %s\n" % way.tags.get("name", "n/a"))
                    f.write("  Highway: %s\n" % way.tags.get("highway", "n/a"))
                    f.write("  Nodes:\n")
                    for node in way.nodes:
                        f.write("    Lat: %f, Lon: %f\n" % (node.lat, node.lon))

    #Keys: Road names
    #Values: List of lists of coordinates (nodes)
    # {Road Name : [[segmentCoordinateTuples]]}
    roadsDict = {}
    with st.spinner("Grouping coordinates by road name"):
        for result in results:
            p += 18/len(results)
            pBar.progress(int(p), "Sorting coordinates...")
            for way in result.ways:
                wayName = way.tags.get("name", "n/a")
                if wayName not in roadsDict.keys():
                    roadsDict[wayName] = []

                wayNodes = []
                for node in way.nodes:
                    wayNodes.append((node.lat, node.lon))

                roadsDict[wayName].append(wayNodes)
    st.success("Grouped coordinates by road name", icon='✔')
    
    if debug:
        print("Grouped by road name")
        with open(fName + "SortedByRoadUncombined.txt", "w") as f:
            for road in roadsDict.keys():
                f.write(road + '\n')
                for nodeList in roadsDict[road]:
                    f.write("Nodes: \n")
                    for node in nodeList:
                        f.write("    Lat: %f, Lon: %f\n" % (node[0], node[1]))
                        f.write('\n')

    with st.spinner("Ordering coordinates into road segments"):
        for road in roadsDict.keys():
            p += 33/len(roadsDict)
            pBar.progress(int(p), "Sorting coordinates...")
            roadsDict[road] = combineSegments(roadsDict[road])
    st.success("Combined coordinates into road segments", icon='✔')

    if debug:
        with open(fName + "CombinedSegments.txt", "w") as f:
            f.write("ROADS FOUND: " + str(len(roadsDict.keys())))
            if home:
                f.write("HOME\n")
                f.write("Palomares Road\n")
                f.write("   Start: " + str(roadsDict["Palomares Road"][0][0]) + ", " + str(roadsDict["Palomares Road"][0][1]) + '\n')
                f.write("   End: " + str(roadsDict["Palomares Road"][-1][0]) + ", " + str(roadsDict["Palomares Road"][-1][1]) + '\n')

                f.write("Dublin Boulevard\n")
                f.write("   Start: " + str(roadsDict["Dublin Boulevard"][0][0]) + ", " + str(roadsDict["Dublin Boulevard"][0][1]) + '\n')
                f.write("   End: " + str(roadsDict["Dublin Boulevard"][-1][0]) + ", " + str(roadsDict["Dublin Boulevard"][-1][1]) + '\n')

            for road in roadsDict.keys():
                f.write(road + '\n')
                for node in roadsDict[road]:
                    f.write("   " + str(node[0]) + ", " + str(node[1]) + '\n')

    return roadsDict

def countSegments(filteredRoadsDict):
    """
        Only works with filtered road dictionaries (filtered for min length, straightness, elevation)
    """
    totalSegments = 0
    for road in filteredRoadsDict.keys():
        for segment in filteredRoadsDict[road]:
            totalSegments += 1

    return totalSegments

def filterRoadsMinLength(roadsDict, minLength, maxLength, pBar, p):
    progAddition = 2
    roadsDictMinLen = {}
    with st.spinner("Slicing roads into segments of requested length"):
        for road in roadsDict.keys():
            p += progAddition/len(roadsDict.keys())
            pBar.progress(int(p), "Filtering for length...")
            longEnoughSegments = []
            segment = []
            startPtr = 0
            endPtr = 1

            while startPtr < len(roadsDict[road]) - 1:
                distance = haversine(roadsDict[road][startPtr][0], roadsDict[road][startPtr][1], roadsDict[road][endPtr][0], roadsDict[road][endPtr][1])
                if distance > maxLength:
                    startPtr += 1
                    endPtr = startPtr
                
                elif distance >= minLength:
                    longEnoughSegments.append((startPtr, endPtr))
                
                if endPtr < len(roadsDict[road]) - 1:
                    endPtr += 1
                else:
                    startPtr += 1
                    
            
            if len(longEnoughSegments) > 0:
                roadsDictMinLen[road] = longEnoughSegments
    st.success("Sliced roads into segments of requested length", icon='✔')

    if debug:
        with open(fName + "MinLength.txt", 'w') as f:
            f.write("FILTERED FOR MIN LENGTH: " + str(minLength) + ", MAX LENGTH: " + str(maxLength) + "\n")
            f.write("SEGMENTS OF MIN LENGTH FOUND: " + str(countSegments(roadsDictMinLen)) + "\n\n")
            for road in roadsDictMinLen.keys():
                f.write(road + '\n')
                for i, segment in enumerate(roadsDictMinLen[road]):
                    roadCoords = roadsDict[road]
                    startCoords = roadCoords[segment[0]]
                    endCoords = roadCoords[segment[1]]
                    f.write("Segment: " + str(i) + '\n')
                    f.write("Length: " + str(haversine(startCoords[0], startCoords[1], endCoords[0], endCoords[1])) + ' ft\n')

    return roadsDictMinLen

def is_straight_road(coordinates, tolerance):
    """
    Determine if a road, represented by a list of coordinates,
    is straight or not.

    Args:
        coordinates (list of tuples): List of (x, y) coordinate pairs
        tolerance (float): Tolerance level for the slope difference

    Returns:
        bool: True if the road is straight, False otherwise
    """
    if len(coordinates) < 3:
        return False  # A road must have at least 3 points to be accurately deemed straight

    # Calculate the slope between the first two points
    x1, y1 = coordinates[0]
    x2, y2 = coordinates[1]
    dy = float((y2 - y1))
    dx = float((x2 - x1))

    if dx == 0:
        dx = .0000000001

    slope = dy/dx
    # Check the slope between each subsequent pair of points
    for i in range(2, len(coordinates)):
        x1, y1 = coordinates[i - 1]
        x2, y2 = coordinates[i]
        newDy = float((y2-y1)) 
        newDx = float((x2-x1))

        if newDx == 0:
            newDx = .00000001

        new_slope = newDy / newDx
        # Check if the absolute difference between slopes is within tolerance
        if abs(new_slope - slope) > tolerance:
            return False

    return True

def filterRoadsStraightness(roadsDict, roadSegmentsMinLen, straightTol, pBar, p):
    progAddition = 2
    roadsDictStraightSegs = {}
    with st.spinner("Checking road segments for straightness"):
        for road in roadSegmentsMinLen.keys():
            p += progAddition/len(roadSegmentsMinLen.keys())
            pBar.progress(int(p), "Filtering for straightness...")
            straightSegs = []
            for segment in roadSegmentsMinLen[road]:
                if is_straight_road(roadsDict[road][segment[0]:segment[1]], straightTol):
                    straightSegs.append(segment)

            if len(straightSegs):
                roadsDictStraightSegs[road] = straightSegs
    st.success("Found roads of requested straightness", icon='✔')


    if debug:
        with open(fName + "Straight.txt" , 'w') as f:
            f.write("STRAIGHT ROADS OF MIN LENGTH:\n")
            f.write("SEGMENTS FOUND: " + str(countSegments(roadsDictStraightSegs)) + "\n\n")
            for road in roadsDictStraightSegs.keys():
                f.write(road + '\n\n')
                for i, segment in enumerate(roadsDictStraightSegs[road]):
                    roadCoords = roadsDict[road]
                    startCoords = roadCoords[segment[0]]
                    endCoords = roadCoords[segment[1]]
                    f.write("Segment: " + str(i) + '\n')
                    # f.write("Start: " + str(roadsDict[road][segment[0]][0]) + ", " + str(roadsDict[road][segment[0]][1]) + '\n')
                    # f.write("End: " + str(roadsDict[road][segment[1]][0]) + ", " + str(roadsDict[road][segment[0]][1]) + '\n')
                    for node in roadsDict[road][segment[0]:segment[-1]]:
                        f.write(str(node[0]) + ", " + str(node[1]) + '\n')
                    f.write("Length: " + str(haversine(startCoords[0], startCoords[1], endCoords[0], endCoords[1])) + ' ft\n\n')

    return roadsDictStraightSegs

def filterDuplicates(roadsDict, filteredRoadsDict, minD):
    for road in filteredRoadsDict.keys():
        i = 0
        while i < len(filteredRoadsDict[road]):
            deleted = False
            segment = filteredRoadsDict[road][i]
            coords = roadsDict[road][segment.start:segment.end]
            for segment2 in filteredRoadsDict[road]:
                contained = True
                if segment2 != segment:
                    startDistance = haversine(roadsDict[road][segment.start][0],
                                              roadsDict[road][segment.start][1],
                                              roadsDict[road][segment2.start][0],
                                              roadsDict[road][segment2.start][1])
            
                    if startDistance < minD * 0.20:
                        deleted = True
                        filteredRoadsDict[road].remove(segment)
                        break

                    coords2 = roadsDict[road][segment2.start:segment2.end]
                    for coord in coords:
                        if coord not in coords2:
                            contained = False

                    if contained:
                        filteredRoadsDict[road].remove(segment)
                        break
            if not deleted:
                i += 1

def get_elevation_data(coordinates):
    now = datetime.now()
    if debug:
        print(now.strftime("%H:%M:%S"))
        print("Getting elevation data")
    time.sleep(1.1)
    elevations = []
    url = "https://api.opentopodata.org/v1/ned10m?locations="

    for coord in coordinates:
        url += str(coord[0]) + "," + str(coord[1]) + "|"

    response = requests.get(url)
    data = response.json()

    for result in data["results"]:
        elevation_meters = result["elevation"]
        if metric:
            elevations.append(elevation_meters)
        else:
            elevation_feet = elevation_meters * 3.28084
            elevations.append(elevation_feet)

    return elevations

def is_flat_road(coordinates, tolerance, elevationData):
    if len(coordinates) < 2: #make sure this is actually a road segment
        return (False, [])
    
    elevations = []
    for c in coordinates:
        elevations.append(elevationData[c])

    startElevation = elevationData[coordinates[0]]

    for i in range(0, len(elevations)):
        for j in range(i + 1, len(elevations)):
            d = elevations[i] - elevations[j]
            if abs(d) > tolerance:
                return (False, [])
    
    return (True, elevations)

def max_elevation_delta(elevations):
    maxD = 0
    for i in range(0, len(elevations)):
        for j in range(i + 1, len(elevations)):
            d = elevations[i] - elevations[j]
            if abs(d) > abs(maxD):
                maxD = d

    return maxD       

def mean_elevation(elevations):
    return statistics.mean(elevations)

def filterRoadsFlatness(roadsDict, filteredRoadSegs, flatTol, pBar, p):
    progAddition1 = 32
    progAddition2 = 5
    roadsDictFlatSegs = {}

    coordinateElevations = {}
    for road in filteredRoadSegs.keys():
        for segment in filteredRoadSegs[road]:
            for coord in roadsDict[road][segment[0]:segment[1]]:
                coordinateElevations[coord] = -1000
    
    toFetch = []
    fetchCounter = 0
    for i, point in enumerate(coordinateElevations.keys()):
            p += progAddition1/len(coordinateElevations.keys())
            pBar.progress(int(p), "Fetching elevation data...")
            toFetch.append(point)
            if i % 99 == 0:
                with st.spinner("Fetching elevation data for coordinates " + str(fetchCounter * 100) + "-" + str((fetchCounter * 100 + 99))):
                    elevations = get_elevation_data(toFetch)
                    for j, elevation in enumerate(elevations):
                        coordinateElevations[toFetch[j]] = elevation
                    toFetch = []
                st.success("Retrieved elevation data for coordinates " + str(fetchCounter * 100) + "-" + str((fetchCounter * 100 + 99)), icon='✔')
                fetchCounter += 1
    
    if len(toFetch) > 0:
        with st.spinner("Fetching elevation data for coordinates " + str(fetchCounter * 100) + "-" + str((fetchCounter * 100) + len(toFetch))):
            elevations = get_elevation_data(toFetch)
            for j, elevation in enumerate(elevations):
                coordinateElevations[toFetch[j]] = elevation
        st.success("Retrieved elevation data for coordinates " + str(fetchCounter * 100) + "-" + str((fetchCounter * 100) + len(toFetch)), icon='✔')
        p += progAddition1/len(coordinateElevations.keys())
        pBar.progress(int(p), "Fetching elevation data...")
    with st.spinner("Checking flatness of road segments"):
        for road in filteredRoadSegs.keys():
            p += progAddition2/len(filteredRoadSegs.keys())
            pBar.progress(int(p), "Checking flatness...")
            flatSegs = []
            for segment in filteredRoadSegs[road]:
                flatness = is_flat_road(roadsDict[road][segment[0]:segment[1]], flatTol, coordinateElevations)
                if flatness[0]:
                    newSegment = ThreeDimSegment(segment[0], segment[1], flatness[1])
                    flatSegs.append(newSegment)

            if len(flatSegs):
                roadsDictFlatSegs[road] = flatSegs
    st.success("Found roads of requested flatness", icon='✔')

    time.sleep(1)
    if debug:
        with open(fName + "Flat.txt" , 'w') as f:
            f.write("FLAT/STRAIGHT ROADS OF MIN LENGTH:\n")
            f.write("SEGMENTS FOUND: " + str(countSegments(roadsDictFlatSegs)) + "\n\n")
            for road in roadsDictFlatSegs.keys():
                f.write(road + '\n\n')
                for i, threeDsegment in enumerate(roadsDictFlatSegs[road]):
                    roadCoords = roadsDict[road]
                    startCoords = roadCoords[threeDsegment.start]
                    endCoords = roadCoords[threeDsegment.end]
                    f.write("Segment: " + str(i) + '\n')
                    for j, node in enumerate(roadsDict[road][threeDsegment.start:threeDsegment.end]):
                        f.write(str(node[0]) + ", " + str(node[1]) + ", " + str(threeDsegment.elevation[j]) + " ft" '\n')
                    f.write("Max elevation change: " + str(max_elevation_delta(threeDsegment.elevation)) + "\n")
                    f.write("Length: " + str(haversine(startCoords[0], startCoords[1], endCoords[0], endCoords[1])) + ' ft\n\n')

    return roadsDictFlatSegs

def countSegments(filteredRoadsDict):
    count = 0
    for road in filteredRoadsDict.keys():
        count += len(filteredRoadsDict[road])

    return count

def findRoads(lat:float, lon:float, radius:float, minDistance:int, maxDistance:int, reqStraightness:float, maxElevChange:int, metricTog:bool):
    progBar = st.progress(0, "Working...")
    prog = 0

    global metric
    metric = metricTog

    if metric:
        unitString = " m"

    else:
        unitString = " ft"
    
    with st.status("Calculating bounding box"):
        print("calculating bbox")
        with st.spinner("Calculating bounding box"):
            bbox = calculate_bounding_box(lat, lon, radius)
        st.success("Bounding box calculated", icon='✔')
        prog += 1
        progBar.progress(prog, "Getting coordinates...")
        print(prog)
    with st.status("Getting local road coordinates"):
        print("getting road data")
        roadsDict = getRoads(bbox[0], bbox[1], bbox[2], bbox[3], progBar, prog)
        prog += 56
        progBar.progress(prog, "Filtering for length...")
    with st.status("Filtering for length"):
        print("filtering for length")
        roadSegmentsMinLen = filterRoadsMinLength(roadsDict, minDistance, maxDistance, progBar, prog)
        prog += 2
        progBar.progress(prog, "Filtering for straightness...")
    with st.status("Filtering for straightness"):
        print("filtering for straightness")
        roadsDictStraightSegs = filterRoadsStraightness(roadsDict, roadSegmentsMinLen, reqStraightness, progBar, prog)
        prog += 2
        progBar.progress(prog, "Filtering for flatness...")
    with st.status("Filtering for flatness", expanded=True) as status:
        print("filtering for flatness")
        roadsDictFlatSegs = filterRoadsFlatness(roadsDict, roadsDictStraightSegs, maxElevChange, progBar, prog)
        print(prog)
        status.update(expanded = False)
        prog += 40
        progBar.progress(prog, "Done!")

    befDupFilter = countSegments(roadsDictFlatSegs)
    filterDuplicates(roadsDict, roadsDictFlatSegs, minDistance)
    aftDupFilter = countSegments(roadsDictFlatSegs)

    if debug:
        print("Before: " + str(befDupFilter))
        print("After: " + str(aftDupFilter))
        print("writing json")

    with st.status("Preparing output"):
        with st.spinner("Writing road data"):
            with open("roadOutput.json", 'w') as f:
                outJsonList = []
                outJsonDetailList = []
                segmentCount = 0
                for k, road in enumerate(roadsDictFlatSegs.keys()):
                    roadCoords = roadsDict[road]
                    for j, segment in enumerate(roadsDictFlatSegs[road]):
                        outDict = {}
                        outDetailDict = {}
                        startCoords = roadCoords[segment.start]
                        endCoords = roadCoords[segment.end]
                        maxElevD = round(max_elevation_delta(segment.elevation), 4)
                        segmentLen = round(haversine(startCoords[0], startCoords[1], endCoords[0], endCoords[1]), 2)
                        outDict["name"] = str(str(road) + "\nMax Elevation Change: " + str(maxElevD) + unitString 
                            + "\nSegment ID: " + str(segmentCount)
                            + "\nLength: " + str(segmentLen) + unitString
                            + "\nStart Coordinates: " + str(startCoords[0]) + ", " + str(startCoords[1])
                            + "\nEnd Coordinates: " + str(endCoords[0]) + ", " + str(endCoords[1])
                            + "\nMean Elevation: " + str(round(mean_elevation(segment.elevation), 4)) + unitString)
                        outDetailDict["Seg ID"] = segmentCount
                        outDetailDict["Road Name"] = str(road)
                        outDetailDict["Max Elev 𝚫"] = maxElevD
                        outDetailDict["Max Elev 𝚫 (abs)"] = abs(maxElevD)
                        outDetailDict["Length (" + unitString + ")"] = segmentLen
                        outDetailDict["Start Coordinates"] = str(startCoords[0]) + ", " + str(startCoords[1])
                        outDetailDict["End Coordinates"] = str(endCoords[0]) + ", " + str(endCoords[1])
                        color = "#" + str(hex(random.randint(0xAAAAAA, 0xFEFEFE)))
                        color = color.replace("0x", '')
                        outDict["color"] = str(color)
                        pathList = []
                        for coord in roadCoords[segment.start:segment.end]:
                            pathList.append([float(coord[1]), float(coord[0])])
                        outDict["path"] = pathList
                        outJson = json.dumps(outDict, indent=4)
                        outJsonDetail = json.dumps(outDetailDict, indent = 4)
                        outJsonList.append(outJson)
                        outJsonDetailList.append(outJsonDetail)
                        segmentCount += 1

                f.write("[\n")
                outString = "["
                outDetailString = "["
                for i, item in enumerate(outJsonList):
                    f.write(item)
                    outString += item
                    outDetailString += outJsonDetailList[i]
                    if (i != len(outJsonList) - 1):
                        f.write(",\n")
                        outString += ","
                        outDetailString += ","                       

                f.write("\n]")
                outString += "]"
                outDetailString += "]"
                        
        st.success("Wrote road data")
        st.session_state.roadData = StringIO(outString)
        st.session_state.roadDetails = StringIO(outDetailString)

    print("done")


if __name__ == "__main__":
    if debug:
        if outlets: 
            originLat = 37.697404
            originLon = -121.848223
        if home:
            originLat = 37.7067134
            originLon = -121.9600764
        if willow:
            originLat = 37.6926971
            originLon = -121.8983176
        if interchange:
            originLat = 37.7009640
            originLon = -121.9226350
    findRoads(originLat, originLon, 2, 1320, 1500, 0.01, 10, False)