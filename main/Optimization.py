# import required packages
import orekit
from orekit.pyhelpers import setup_orekit_curdir, download_orekit_data_curdir

from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint
from org.orekit.time import TimeScalesFactory, AbsoluteDate, DateComponents, TimeComponents
from org.orekit.utils import Constants, IERSConventions, PVCoordinates
from org.orekit.propagation.analytical.tle import TLE,TLEPropagator
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from math import radians, pi, degrees
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from org.orekit.propagation.analytical import EcksteinHechlerPropagator, KeplerianPropagator

from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.geometry.fov import CircularFieldOfView, FieldOfView
from org.orekit.propagation.events import FieldOfViewDetector, ElevationDetector, BooleanDetector, EventsLogger
from datetime import datetime, timedelta
from org.orekit.propagation.events.handlers import ContinueOnEvent

from org.orekit.attitudes import NadirPointing

# from math import degrees
from org.orekit.frames import Transform
from org.hipparchus.geometry import Vector

import numpy as np

import json
from shapely.geometry import Point, Polygon, GeometryCollection
from shapely.geometry import MultiPoint, MultiPolygon
from shapely.geometry import mapping
from shapely.geometry import shape
from pyproj import Geod
from shapely import wkt

from java.util import List
from datetime import datetime, timedelta
import os

from deap import base, creator, tools, algorithms


# PROJECT_DIR = os.path.join(os.path.sep, 'projects', os.environ['USER'], 'my_project')
# DATA_IN = os.path.join(PROJECT_DIR, 'data_in')
# DATA_OUT = os.path.join(PROJECT_DIR, 'data_out', os.environ['SLURM_JOB_ID'])



# import defined functions
# from tleParser import parse_tle 
# from SatelliteState import get_satelliteStates, propagator_to_statePosition

# default tle for this project
# IRIDIUM 178
line1="1 56729U 23068V   24154.47728178  .00000642  00000-0  81118-4 0  9993"
line2= "2 56729  86.5773 339.2187 0002054  81.2330 278.9124 14.80208184 56200"


westernCentral_partArea = [[
            147.17515827251634,
            -15.640078251023056
          ],
          [
            147.94567963567386,
            -12.187583366380892
          ],
          [
            152.47710504377255,
            -11.812564622640394
          ],
          [
            157.41096167261662,
            -12.652514202747142
          ],
          [
            160.04819327494533,
            -15.759696732458849
          ],
          [
            158.98519773901563,
            -19.70058922777251
          ],
          [
            154.059202587223,
            -19.649308459193833
          ],
          [
            150.1129330003555,
            -18.105241744613366
          ],
          [
            147.2062653984351,
            -15.552334895009977
          ]]


#  function to parse the tle to keplarian elements
def parse_tle(line1, line2): 
    #this function will take the tle information (line1, line2) and return the dict value of keplerian orbital parameters formated as
    """
       {
       'mean_motion': ,
       'eccentricity': ,
       'inclination': ,
       'right_ascension': ,
       'argument_of_perigee': ,
       'mean_anomaly': ,
       'semi_major_axis': 
       }
    """

    keplerian = {}
    
    keplerian['mean_motion'] = float(line2[52:63])
    keplerian['eccentricity'] = float('0.' + line2[26:33])
    keplerian['inclination'] = float(line2[8:16])
    keplerian['right_ascension'] = float(line2[17:25])
    keplerian['argument_of_perigee'] = float(line2[34:42])
    keplerian['mean_anomaly'] = float(line2[43:51])
        # Calculate semi-major axis (a) using mean motion (n) and Earth's gravitational constant (mu)
    mu = 398600.4418  # Earth's gravitational constant in km^3/s^2      which is approx 3.986004418(8)×10^14 m^3⋅s^−2
    n = keplerian['mean_motion'] * 2 * pi / 86400  # Convert mean motion to radians/minute

    # The semi-major axis of an orbit can be found using Kepler's Third Law, which states that the square of the orbital period (T) of a planet is directly proportional to the cube of the semi-major axis (a) of its orbit. 
    # This can be written as: T^2 = k * a^3 where k is a constant of proportionality
    a = (mu / (n ** 2)) ** (1 / 3)  # Calculate semi-major axis in kilometers
    keplerian['semi_major_axis'] = a
    
    return keplerian


"""
THIS get_satelliteStates function is to get the [SpacecraftState<>, SpacecraftState<>,......] ie, the state of satellite with time/epoch, 
"""
def get_satelliteStates(propagator, shiftedDuration_sec, propagator_duration):
    # this will take propagator and return list of "SpacecraftState"
    # irridium original tle orbit
    currentStateIrr = propagator.getInitialState()
    initialDateIrr = currentStateIrr.getDate()   #as this this return the AbsoluteDate

    endDateIrr = initialDateIrr.shiftedBy(propagator_duration)
    
    # print("initial date of TLE orbit irridium",initialDateIrr)
    
    statesIrri = []
    
    while (currentStateIrr.getDate().compareTo(endDateIrr)<=0):
        statesIrri.append(currentStateIrr)
        currentStateIrr = propagator.propagate(currentStateIrr.getDate().shiftedBy(shiftedDuration_sec))  #shifted every 60 sec
        # print(currentStateIrr.getDate())
    return statesIrri


"""
For generation of czml file, we need the satellite position information in format of [clock, px, py, pz] and that should also be for multiple satellite or propagator, so,
we have this function to get the position, 
"""

# expected return values should be as [ ([[c, px,py,pz],[c1,px1,py1,pz1]], start_date, end_date), ([[]],ds,de), ....]
# which take list of propagator, and gives the position at any time, for each propagator/satellte
def propagator_to_statePosition(listOfPropagator, shiftedDuration_sec, propagato_duration):

    listOfSatPosition = []
    for id,propagator_i in enumerate(listOfPropagator):
        # propagator = propagator
        listOfStates = []
        # propagator_0 = listOfPropagator[0]
        currentState = propagator_i.getInitialState()
        start_date = currentState.getDate()
        end_date = start_date.shiftedBy(propagato_duration)
        
        satStates = get_satelliteStates(propagator_i, shiftedDuration_sec, propagato_duration)
        listOfStates.append((satStates,start_date, end_date))
        # listOfStates.append(get_satelliteStates(listOfPropagator[1], shiftDuration_sec))
        
        for state_id, (satState,start_date,end_date) in enumerate(listOfStates):
            # Define satellite points and trajectory
            satellite_id = state_id
            startDate = start_date
            endDate = end_date
            duration = 0.0
            trajectory = []
            for i,state in enumerate(satState):
                position = state.getPVCoordinates().getPosition().toArray() # extract the position from Vector3D into array
                trajectory.append((duration,position[0], position[1],position[2]))
                duration = duration+shiftedDuration_sec # this duration shift is for clock count, as we have a propagator shifted by same shiftedDuration_sec
                
            listOfSatPosition.append((trajectory,start_date,end_date))
    # print (datetime.now() - startTime)
    return listOfSatPosition



def create_orbit_propagator(semi_axis_km, ecc, i_degree, a_p_degree, raan_degree, ma_degree, epoch):
    orbit = KeplerianOrbit(semi_axis_km * 1000, ecc, radians(i_degree), radians(a_p_degree), radians(raan_degree), 
                           radians(ma_degree), PositionAngleType.MEAN, FramesFactory.getEME2000(), epoch, Constants.WGS84_EARTH_MU)
    return KeplerianPropagator(orbit)

def get_footprint(propagator, date, earth):
    spaceCraft_state = propagator.propagate(date)
    
    # satellite position 
    statePoint = earth.transform(spaceCraft_state.getPosition(), spaceCraft_state.getFrame(), spaceCraft_state.getDate())   # this will transform orbit frame TEME to ITRF earthFrame for point
    satlon = degrees(statePoint.getLongitude())
    satlat =degrees(statePoint.getLatitude())
    satAlt =statePoint.getAltitude()
    
    #Caculate the los in ECF frame
    target_geodeticPoint = GeodeticPoint(radians(satlat), radians(satlon), 0.0);
    targetEarth = earth.transform(target_geodeticPoint)
    los_sat2target_ecf   = targetEarth.subtract(spaceCraft_state.getPVCoordinates(earthFrame).getPosition())
    
    #get the transform from J2000 to ECF
    inertToBody = spaceCraft_state.getFrame().getTransformTo(earth.getBodyFrame(), date)
    
    #get the transform from satellite body frame to ECF frame
    satBody2ecf = Transform(date, spaceCraft_state.toTransform().getInverse(), inertToBody)
    #transform the los into the satellite body frame
    losPointing_body=satBody2ecf.getInverse().getRotation().applyTo(Vector.cast_(los_sat2target_ecf).normalize())
    
    
    fov = CircularFieldOfView(losPointing_body, radians(25.0),  radians(0.0))       #just taking the center as nadir point as it is giving a  circular footprint, need to look further about the axes of spacecraft and sensor 
    footPrintList=fov.getFootprint(satBody2ecf, earth, 0.1)
    return footPrintList, satlon, satlat

def casting_footprintList(footPrintList):
    footprint = list(footPrintList)
    footprint = [list(List.cast_(points)) for points in footprint]
    listOfPoints = []
    for points in footprint:
            for point in points:
                    
                    point = GeodeticPoint.cast_(point)
                    listOfPoints.append(point)
    return listOfPoints


def detectAreaOfInterest(listOfPoints, areaOfInterest, date):
    # we will be using the exciting function of polygon here to find the interesection
    polygon_areaOfInterest = Polygon(areaOfInterest)
    polygon_footPrint = Polygon([[degrees(p.getLongitude()), degrees(p.getLatitude())] for p in listOfPoints])
    #  to avoid the TopologyException: side location conflict, may be due to a bowtie, in polygon getting this issue, it is recommended to solve this exception using shapely buffer(0) method as
    polygon_footPrint = polygon_footPrint.buffer(0)
    # print("area of footprint simply calculated: ", polygon_footPrint.area)

    if(polygon_footPrint.area > 450):
        """
        On checking manually, when the satellite is at the edge of earth 2d plane, circular polygon shape got elongated as cylindrical,
        as spherical shape of earth surface is unfold to 2d plane, this issue was giving footprint area for 2d plane surface too much, so needed to filter as exception.
        
        """
        return (False, 0.0, 0.0)
    # Calculate the intersection (overlap) of the polygons
    intersect = polygon_areaOfInterest.intersection(polygon_footPrint)
    # print(intersect)
    area_covered_simply = intersect.area
    # print("percentage of area covered by footprint on area of interset: ",area_covered_simply)

    # we can also calculate the overlap/intersect area as the actual earth surface area in sq meter as
    interset_str = str(intersect)
    interset_poly = wkt.loads(interset_str)
    geod = Geod(ellps="WGS84")
    geod_area_sq_m = abs(geod.geometry_area_perimeter(interset_poly)[0])
    
    # print(f"itersect Geodesic area: {geod_area_sq_m:.3f} m^2")

        # print( polygon_areaOfInterest.intersects(polygon_footPrint))
    if(polygon_areaOfInterest.intersects(polygon_footPrint) == True):
        # print( polygon_areaOfInterest.intersects(polygon_footPrint))
        return (True, area_covered_simply, geod_area_sq_m)
    else:
        return (False, 0.0, 0.0)


# taking propagator, period of time, and detectAreaofInterest together with small steps to calculate the dwell period of satellite on area of interest

def loopWholeFunction(propagator, startDate, duration,areaOfInterest, earth):
    # codeblockStartTime = datetime.now()

    setFlag = 0
    notedDate = []
    calculatedDwellTime = []
    dwellTimeInterval = []
    date = startDate
    while (date.compareTo(startDate.shiftedBy(duration)))<=0:
        # print(date)
        footprintlist, lon_craft, lat_craft = get_footprint(propagator, date, earth)
        # print("craft position ", lon_craft, lat_craft)
        castedFootPrintlist = casting_footprintList(footprintlist)
        detection, area, geod_area_sq_meter = detectAreaOfInterest(castedFootPrintlist,areaOfInterest, date)
        if (detection == True):
            # print("detected")
            # print(f"percentage of area detected by sensor footPrint: {geod_area_sq_meter}")
            # print(date)
            notedDate.append(date)
            date = date.shiftedBy(10.0)
            setFlag = 1
        else:   
                
            if (len(notedDate)>0 and (setFlag == 1)):
                date = date.shiftedBy(60.0*40.0) # this will skip 40 min for generating footprint and detection, cause it is sure for that, there won't be any detection early that early
                # just by skiping time we get the benefit of more than 1 sec computation cost
                # print("time shifted by 40 minutes")
                setFlag = 0
                
                startDetection = notedDate[0]
                endDetection = notedDate[-1]
                dwellTimeInterval.append([startDetection.toString(utc),endDetection.toString(utc)])
                calculatedDwellTime.append(endDetection.durationFrom(startDetection))
            # print(date)    
            date = date.shiftedBy(20.0)
            notedDate.clear()
            
    # calculatedDwelltime =  find the time difference first and last time of notedDate
    # print(datetime.now() - codeblockStartTime)
    return calculatedDwellTime, dwellTimeInterval
    
   
# function to calculate the revisit time using the dwell time interval data

#  now we don't have to deal with the Absolutedate or 'Z' anymore, cause I already stored the time information as utc string 
def parse_datetime(date_str):
    return datetime.fromisoformat(date_str.replace('Z', '+00:00'))


def calculate_revisit_time(dwell_timeInterval_data):
    cummulative_revisit_time = timedelta(seconds=0)
    # Sort the data based on the first element of the inner list
    sorted_data = sorted(dwell_timeInterval_data, key=lambda x: (x[1][0]))

    # Calculate time differences between consecutive intervals
    time_differences = []
    for i in range(len(sorted_data) - 1):
        end_current = parse_datetime(sorted_data[i][1][1])
        start_next = parse_datetime(sorted_data[i + 1][1][0])
        
        difference = start_next - end_current
        cummulative_revisit_time = cummulative_revisit_time + difference
        time_differences.append((sorted_data[i], sorted_data[i + 1], difference))
    
    # when there is empty list of dwell time interval data, or ie, there is no detection, we are getting cummulative dwell time as 0 sec, which is the optimal minimum so, resolving this error 
    if (cummulative_revisit_time == timedelta(seconds = 0)):
        cummulative_revisit_time = timedelta(seconds = 865000)
        return time_differences, cummulative_revisit_time
    else:
        return time_differences, cummulative_revisit_time


def calculate_cumm_revisit_dwell_time(satellites, target_location, simulation_duration, epoch, earth):
    
    dwell_time_sats = []
    dwell_time_intervals = []
    revisit_time_sats = []
    
    for id, propagator in enumerate(satellites):
        
        calDT, dtInt = loopWholeFunction(propagator,epoch, simulation_duration, target_location, earth)
        
        dwell_time_sats.append((id,calDT))
        dwell_time_intervals.append((id,dtInt))
    
    # flatten the list of dwell time interval list data first
    desired_dtIntOfSat = [(index, value) for index, values in dwell_time_intervals for value in values]
    # calcualte revisit time, and cummulative revisit time
    revist_time_differences, cumm_revisit_time = calculate_revisit_time(desired_dtIntOfSat)
   
    
    # this will gives the cummulative dwell time
    total_sum = sum(sum(values) for _, values in dwell_time_sats)
    cummulative_dwell_time = timedelta(seconds = total_sum)
   
    return cumm_revisit_time,cummulative_dwell_time




def main():

    # from org.orekit.sensors import FieldOfView, ImagingSensor, ObservationSensor
    #  /// initialization of package


    mytle = TLE(line1, line2)
    kepOrbEle = parse_tle(line1,line2)

    print(kepOrbEle)

    epoch = AbsoluteDate(2024, 1, 1, 0, 0, 0.000, TimeScalesFactory.getUTC())
    fitness_values = []

    def create_walker_delta_constellation(num_planes, num_sats_per_plane, inclinations_deg):
        satellites = []     #for list of propagator
        mu = Constants.WGS84_EARTH_MU
        # inclination = radians(inclination_deg)

        for plane in range(num_planes):
            inclination =inclinations_deg[plane]    #degrees
            for sat in range(num_sats_per_plane):
                raan = degrees(2 * np.pi * plane / num_planes)   #degrees
                mean_anomaly = degrees(2 * np.pi * sat / num_sats_per_plane)    #degrees
                propagator = create_orbit_propagator(kepOrbEle.get("semi_major_axis"), 0.001, inclination,
                                            kepOrbEle.get("argument_of_perigee"),raan,
                                            mean_anomaly,
                                            epoch)
                satellites.append(propagator)
        
        return satellites      
    
    # Define the evaluation function for the GA
    def evaluate_revisit_time(individual):
        num_planes = len(individual)
        num_sats_per_plane = 4  # Example number of satellites per plane

        satellites = create_walker_delta_constellation(num_planes, num_sats_per_plane, individual)
        target_location = westernCentral_partArea
        simulation_duration = 86400.0  # 1 day in seconds

        codeblockTime= datetime.now()

        cumm_revisit_time, cumm_dwell_time = calculate_cumm_revisit_dwell_time(satellites, target_location, simulation_duration,epoch, earth)

        print(f"total time it takes to cal cumm re time {(datetime.now() - codeblockTime)}")
        print("cummulative revisit time: ",cumm_revisit_time)

        fitness_values.append(cumm_revisit_time)
        
        return cumm_revisit_time,



    def evaluate_dwell_time(individual):
        num_planes = len(individual)
        num_sats_per_plane = 4  # Example number of satellites per plane

        satellites = create_walker_delta_constellation(num_planes, num_sats_per_plane, individual)
        target_location = westernCentral_partArea
        simulation_duration = 86400.0  # 1 day in seconds

        codeblockTime= datetime.now()

        cumm_revisit_time, cumm_dwell_time = calculate_cumm_revisit_dwell_time(satellites, target_location, simulation_duration, epoch, earth)

        print(f"total time it takes to cal cumm re time {(datetime.now() - codeblockTime)}")
        # print("cummulative revisit time: ",cumm_revisit_time)
        return cumm_dwell_time,



    # Create the fitness and individual classes
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()

    # Attribute generator: Define inclination angle for each orbital plane in degrees
    def random_inclination():
        return np.random.uniform(0, 180)
    # Structure initializers: Each individual is a list of inclinations for each plane
    NUM_PLANES = 3  # Example number of orbital planes




    toolbox.register("individual", tools.initRepeat, creator.Individual, random_inclination, n=NUM_PLANES)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("evaluate", evaluate_revisit_time)
    toolbox.register("mate", tools.cxBlend, alpha=0.5)
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1.0, indpb=0.2)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # Run the Genetic Algorithm
    population = toolbox.population(n=2)
    algorithms.eaSimple(population, toolbox, cxpb=0.9, mutpb=0.1, ngen=2, stats=None, halloffame=None, verbose=True)

    # Convert timedelta to hours
    hours = [td.total_seconds() / 3600 for td in fitness_values]
    plt.bar(range(len(hours)), hours, color='blue')
    plt.xlabel('Event Index')
    plt.ylabel('Time (Hours)')
    plt.title('Time Intervals in Hours')
    plt.xticks(range(len(hours)), rotation=90)
    plt.grid(True)
    # Show bar plot
    plt.tight_layout()
    plt.savefig(os.path.join('ConstellationProj/masterResearchProj/main/DATA_OUT', 'timeInterval'))


    # Extract and print the best solution
    best_individual = tools.selBest(population, k=1)[0]
    best_fitness = best_individual.fitness.values[0]
    print(f"Best individual: {best_individual} and best fitness values: {best_fitness}")

# hook
if __name__== '__main__':
    vm  = orekit.initVM()
    download_orekit_data_curdir()    # download the orkeit-data for once 
    setup_orekit_curdir()
    print("ok")
    utc = TimeScalesFactory.getUTC()
    # # //set inertial frame (ECI)
    inertialFrame = FramesFactory.getEME2000()
    # non-inertial Frame
    earthFrame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    # earth model
    earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING, earthFrame)

    main()