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
import numpy as np
from deap.tools import sortNondominated
from deap.tools.emo import assignCrowdingDist



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

    def create_walker_delta_constellation_Multi(num_planes, num_sats_per_plane, inclinations_deg, raan_deg):
        satellites = []     #for list of propagator
        mu = Constants.WGS84_EARTH_MU
        
        # inclination = radians(inclination_deg)

        for plane in range(num_planes):
            inclination =inclinations_deg[plane]    #degrees
            raan = raan_deg[plane]
            for sat in range(num_sats_per_plane):
                # raan = degrees(2 * np.pi * plane / num_planes)   #degrees
                mean_anomaly = degrees(2 * np.pi * sat / num_sats_per_plane)    #degrees
                propagator = create_orbit_propagator(kepOrbEle.get("semi_major_axis"), 0.001, inclination,
                                            kepOrbEle.get("argument_of_perigee"),raan,
                                            mean_anomaly,
                                            epoch)
                satellites.append(propagator)
        return satellites      
    
    def evaluate_both(individual):

        num_planes = len(individual[0])
        num_sats_per_plane = 4  # Example number of satellites per plane
        print("individual:  ", individual)
        inclinations_list = individual[0]    #[plane[0] for plane in individual]  # Extract inclinations
        raans_list = individual[1]    #[plane[1] for plane in individual]  # Extract RAANs
        
        # print("inclination list : ", inclinations_list)
        
        satellites = create_walker_delta_constellation_Multi(num_planes, num_sats_per_plane, inclinations_list, raans_list)
        target_location = westernCentral_partArea
        simulation_duration = 86400.0  # 1 day in seconds

        codeblockTime= datetime.now()

        cumm_revisit_time, cumm_dwell_time = calculate_cumm_revisit_dwell_time(satellites, target_location, simulation_duration,epoch,earth)

        print(f"total time it takes to cal cumm re time {(datetime.now() - codeblockTime)}")
        # print("cummulative revisit time: ",cumm_revisit_time)
        fitness_values.append((cumm_revisit_time, cumm_dwell_time))
        cumm_revisit_time_seconds = cumm_revisit_time.total_seconds()
        cumm_dwell_time_seconds = cumm_dwell_time.total_seconds()
        return cumm_revisit_time_seconds, - cumm_dwell_time_seconds     # set the dwell time -ve as nsga3 minimize the objectives, so setting negative will maximize it


    def generate_reference_points(nobj, p):
        """
        Generate reference points for NSGA-III based on the number of objectives (nobj) and p divisions.
        """
        def recursive_combinations(n, left, current_comb, combinations):
            if n == 0:
                combinations.append(current_comb + [left])
            else:
                for i in range(left + 1):
                    recursive_combinations(n - 1, left - i, current_comb + [i], combinations)

        combinations = []
        recursive_combinations(nobj - 1, p, [], combinations)
        combinations = np.array(combinations) / p
        return combinations

    # Generate reference points for 2 objectives with 12 divisions
    reference_points = generate_reference_points(nobj=2, p=12)

    # Plotting reference points
    def plot_reference_points(reference_points):
        plt.scatter(reference_points[:, 0], reference_points[:, 1], color='green', marker='x', label='Reference Points')
        plt.savefig(os.path.join('DATA_OUT', 'reference_points_nsga3_simple'))
    # Call this after generating reference points
    plot_reference_points(reference_points)




    def select_nsga3(population, k, reference_points):
        """
        Custom NSGA-III selection process with reference points.
        """
        # Sort the population into different fronts
        pareto_fronts = sortNondominated(population, k, first_front_only=False)

        # Assign crowding distance to individuals in the first Pareto front
        for front in pareto_fronts:
            assignCrowdingDist(front)

        # Selection step: Pick the best individuals according to NSGA-III
        selected = []
        for front in pareto_fronts:
            selected.extend(front)
            if len(selected) >= k:
                break

        return selected[:k]
    
    # Extract fitness values from the final population
    def plot_population(population):
        revisit_times = [ind.fitness.values[0] for ind in population]  # Revisit time (minimized)
        print(revisit_times)
        dwell_times = [ind.fitness.values[1] for ind in population]    # Dwell time (maximized, -ve values)
        print(dwell_times)
        # Plot population fitness
        plt.figure(figsize=(10, 6))
        plt.scatter(revisit_times, [-d for d in dwell_times], color='blue', label='Final Population')
        plt.title('Final Population Fitness')
        plt.xlabel('Revisit Time (seconds)')
        plt.ylabel('Dwell Time (seconds)')
        plt.grid(True)
        # plt.legend()
        # plt.show()
        plt.savefig(os.path.join('DATA_OUT', 'final_population_fitness_nsga3_simple'))

    def plot_fitness_values_per_individual(fitness_values):
        # Plotting the fitness values of cum dt and rt, with each individual 
        # fitness_values has data as [(rt1, dt1), (rt2,dt2), (rt3,dt3).....] for each individual
        # Convert timedelta to hours
        revisit_time_hours = [td[0].total_seconds() / 3600 for td in fitness_values]
        dwell_time_hours = [td[1].total_seconds() / 3600 for td in fitness_values]

        # Generate a line plot for comparison
        plt.figure(figsize=(14, 8))
        plt.plot(range(len(revisit_time_hours)), revisit_time_hours, marker='o', linestyle='-', color='red', label="revisit Time Intervals")
        plt.xlabel('Event Index')
        plt.ylabel('Time (Hours)')
        plt.title('revisit Time Intervals over Events')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join('DATA_OUT', 'revisit_time_per_individuals_nsga3_simple'))

        # Generate a line plot for comparison
        plt.figure(figsize=(14, 8))
        plt.plot(range(len(dwell_time_hours)), dwell_time_hours, marker='o', linestyle='-', color='red', label="dwell Time Intervals")
        plt.xlabel('Event Index')
        plt.ylabel('Time (Hours)')
        plt.title('dwell Time Intervals over Events')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join('DATA_OUT', 'dwell_time_per_individuals_nsga3_simple'))

        # plt.show()


    # Define a multi-objective fitness function
    creator.create("FitnessMulti", base.Fitness, weights=(-1.0, -1.0))     #minimizing the both objectives

    # Define an individual containing lists of inclinations and RAANs
    creator.create("Individual", list, fitness=creator.FitnessMulti)

    toolbox = base.Toolbox()

    # Random inclination and RAAN for each plane
    def random_inclination():
        return np.random.uniform(0, 180)

    def random_raan():
        return np.random.uniform(0, 360)

    # Individual representation: Inclination and RAAN for each plane
    def create_individual(num_planes):
        """
        Create an individual where each individual consists of both inclination and RAAN for each plane.
        """
        inclination_list = [random_inclination() for _ in range(num_planes)]
        raan_list = [random_raan() for _ in range(num_planes)]
        return creator.Individual([inclination_list, raan_list])

    # Number of planes in the constellation
    NUM_PLANES = 3

    toolbox.register("individual", create_individual, NUM_PLANES)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)


    # Corrected crossover function (element-wise blending)
    def cxBlendPlanes(ind1, ind2, alpha=0.5):
        print("ind1 : ",ind1)
        print("ind2 : ", ind2)
        for i in range(len(ind1[0])):  # Iterate over each plane
            # Blending inclination
            ind1[0][i] = (1 - alpha) * ind1[0][i] + alpha * ind2[0][i]
            ind2[0][i] = alpha * ind1[0][i] + (1 - alpha) * ind2[0][i]
            
            # Blending RAAN
            ind1[1][i] = (1 - alpha) * ind1[1][i] + alpha * ind2[1][i]
            ind2[1][i] = alpha * ind1[1][i] + (1 - alpha) * ind2[1][i]
        
        return ind1, ind2

    # Custom mutation function (element-wise mutation)
    def mutGaussianPlanes(individual, mu, sigma, indpb):
        for i in range(len(individual[0])):  # Iterate over each plane
            if np.random.random() < indpb:
                inclination = individual[0][i]
                raan = individual[1][i]
                
                # Mutate inclination and RAAN
                mutated_inclination = np.clip(inclination + np.random.normal(mu, sigma), 0, 180)
                mutated_raan = np.clip(raan + np.random.normal(mu, sigma), 0, 360)
                
                # Update the individual with the mutated values
                individual[0][i] = mutated_inclination
                individual[1][i] = mutated_raan
        print("mutation individual: ", individual)
        return individual,

    # Dynamic crossover and mutation adjustment based on generation
    def dynamic_crossover_mutation(gen, ngen, Pc_start=0.7, Pm_start=0.2, Pm_end=0.4):
        Pc = Pc_start * (1 - gen / ngen)  # Crossover rate decreases over generations
        Pm = Pm_start + (Pm_end - Pm_start) * (gen / ngen)  # Mutation rate increases over generations
        return Pc, Pm


    toolbox.register("evaluate", evaluate_both)

    # Register genetic operators
    toolbox.register("mate", cxBlendPlanes)   # Custom crossover function
    toolbox.register("mutate", mutGaussianPlanes, mu=0, sigma=1.0, indpb=0.2)  # Custom mutation function

    # Register NSGA-III selection
    reference_points = generate_reference_points(nobj=2, p=12)
    toolbox.register("select", select_nsga3, reference_points=reference_points)

    # Genetic algorithm with dynamic crossover and mutation probabilities
    def run_ga_dynamic(pop_size=20, ngen=10):
        population = toolbox.population(n=pop_size)
        
        for gen in range(ngen):
            Pc, Pm = dynamic_crossover_mutation(gen, ngen)  # Dynamic crossover and mutation probabilities
            
            offspring = toolbox.select(population, len(population))
            offspring = list(map(toolbox.clone, offspring))
            
            # Apply crossover (with dynamic Pc)
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if np.random.random() < Pc:
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values
            
            # Apply mutation (with dynamic Pm)
            for mutant in offspring:
                if np.random.random() < Pm:
                    toolbox.mutate(mutant)
                    del mutant.fitness.values
            
            # Evaluate individuals with invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            
            # Replace the population with the new offspring
            population[:] = offspring
            
            # Print best individual of the current generation
            best_ind = tools.selBest(population, 1)[0]
            print(f"Generation {gen}: Best Individual = {best_ind}, Fitness = {best_ind.fitness.values}")
            
        return population
        
    nsga_pop = run_ga_dynamic(pop_size=20, ngen=2)  # Example population of 20 and 10 generations

    # # Create population
    # population = toolbox.population(n=4)

    # # Run NSGA-III (similar to NSGA-II)
    # algorithms.eaMuPlusLambda(population, toolbox, mu=10, lambda_=20, cxpb=0.7, mutpb=0.2, ngen=2, stats=None, halloffame=None, verbose=True)


    ####################################################
    # plot the final population 
    plot_population(nsga_pop)
    # plot the fitness_values per individuals event
    plot_fitness_values_per_individual(fitness_values)

    # Extract the Pareto front
    pareto_front = tools.sortNondominated(nsga_pop, len(nsga_pop), first_front_only=True)[0]

    # Plot the Pareto front
    revisit_times = [ind.fitness.values[0] for ind in pareto_front]
    dwell_times = [ind.fitness.values[1] for ind in pareto_front]

    plt.figure(figsize=(10, 6))
    plt.scatter(revisit_times, dwell_times, color='blue', label="Pareto Front")
    plt.title('Pareto Front of NSGA-III')
    plt.xlabel('Revisit Time')
    plt.ylabel('Dwell Time (seconds)')
    plt.grid(True)
    # plt.legend()
    # plt.show()
    plt.savefig(os.path.join('DATA_OUT', 'pareto_front_nsga3_simple'))


    # Get the top 10 best individuals from the Pareto front based on fitness values
    top_10_best = sorted(pareto_front, key=lambda ind: ind.fitness.values)[:10]
    # Print the top 10 individuals with their fitness values (revisit time and dwell time)
    top_10_fitness_values = [(ind.fitness.values[0], ind.fitness.values[1]) for ind in top_10_best]

    # Display the top 10 best fitness values
    for i, fitness in enumerate(top_10_fitness_values, 1):
        print(f"Individual {i}: Revisit Time = {fitness[0]}, Dwell Time = {fitness[1]} seconds")


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