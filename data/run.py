import math
import argparse
import gzip

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--graph-csv', '-g', type=str, help='Path to the graph CSV file')
parser.add_argument('--check-file', '-c', type=str, help='Path to the check file')
parser.add_argument('--endpoints', '-e', type=str, help='Path to the endpoints file')
parser.add_argument('--output', '-o', type=str, help='Path to the output Slack file')
parser.add_argument('--nsigma', '-n', type=int, default=3, help='Value for N sigma (default: 3)')
parser.add_argument('--ignore-sigma', '-i', action='store_true', help='Ignore sigma')
parser.add_argument('--debug', '-d', type=int, default=0, help='Debug mode (default: 0)')

# Parse command-line arguments
args = parser.parse_args()

# Assign parsed arguments to variables
graphCSVFile = args.graph_csv
checkFile = args.check_file
endpointsFile = args.endpoints
outputSlackFile = args.output
N = args.nsigma
ignoreSigma = args.ignore_sigma
debug = args.debug

# Set clock value
clock = 10



def open_file(filename, mode):
    if mode == 'read' or mode == 'r':
        if filename.endswith(".gz"):
            return gzip.open(filename,"rt")
        else:
            return open(filename,"wt")
    elif mode == "write" or mode == "w":
        if filename.endswith(".gz"):
            return gzip.open(filename,"wt")
        else:
            return open(filename,"w")
    else:
        raise ValueError(f"Unsupported mode {mode} in open_file\n")


delay = {}

def add_delay_arc(line):
    global ignoreSigma
    from_node, to_node, sense, max_rise_mean, max_rise_sigma, max_fall_mean, max_fall_sigma, min_rise_mean, min_rise_sigma, min_fall_mean, min_fall_sigma = map(str.strip, line.split(','))
    key = f"{from_node}:{to_node}"

    if ignoreSigma:
        max_rise_sigma = 0
        max_fall_sigma = 0
        min_rise_sigma = 0
        min_fall_sigma = 0

    arc = {
        'sense': sense,
        'mean': [max_rise_mean, max_fall_mean, min_rise_mean, min_fall_mean],
        'sigma': [max_rise_sigma, max_fall_sigma, min_rise_sigma, min_fall_sigma]
    }

    delay.setdefault(key, []).append(arc)

def get_delay_arc(from_node, to_node):
    key = f"{from_node}:{to_node}"
    return delay.get(key, [])


def is_valid_arrival(arr):
    if (
        arr['mean'][0] != -1e99
        or arr['mean'][1] != -1e99
        or arr['mean'][2] != 1e99
        or arr['mean'][3] != 1e99
    ):
        return 1
    else:
        return 0

def create_arrival():
    arr = {'pos_edge': 1}
    mean = [-1e99, -1e99, 1e99, 1e99]
    sigma = [0, 0, 0, 0]
    arr['mean'] = mean
    arr['sigma'] = sigma
    return arr

def init_arrival(arr, pos_edge):
    arr['pos_edge'] = pos_edge
    if pos_edge:
        arr['mean'][0] = 0
        arr['mean'][2] = 0
    else:
        arr['mean'][1] = 0
        arr['mean'][3] = 0

def init_root_arrival(roots):
    for pin in roots:
        arrivals = []
        pos_arr = create_arrival()
        init_arrival(pos_arr, True)
        neg_arr = create_arrival()
        init_arrival(neg_arr, False)
        arrivals.extend([pos_arr, neg_arr])
        set_arrival_times(pin, arrivals)

def find_match_arrival(arr, arrs):
    for a in arrs:
        if arr['pos_edge'] == a['pos_edge']:
            return a
    arrs.append(arr)
    return arr

def add(old_mean, old_sigma, delay_mean, delay_sigma):
    if old_mean != 1e99 and old_mean != -1e99:
        new_mean = old_mean + delay_mean
        new_sigma = math.sqrt(old_sigma**2 + delay_sigma**2)
        return new_mean, new_sigma
    else:
        return old_mean, 0

def add_normal(arr, arr_index, delay, new_arr, new_index):  ## perhaps wrong here 
    max_index_offset = 0

    old_mean = arr['mean'][arr_index + max_index_offset]
    old_sigma = arr['sigma'][arr_index + max_index_offset]
    delay_mean = delay['mean'][new_index + max_index_offset]
    delay_sigma = delay['sigma'][new_index + max_index_offset]
    new_mean, new_sigma = add(old_mean, old_sigma, float(delay_mean), float(delay_sigma))
    new_arr['mean'][new_index + max_index_offset] = new_mean
    new_arr['sigma'][new_index + max_index_offset] = new_sigma

    max_index_offset = 2

    old_mean = arr['mean'][arr_index + max_index_offset]
    old_sigma = arr['sigma'][arr_index + max_index_offset]
    delay_mean = delay['mean'][new_index + max_index_offset]
    delay_sigma = delay['sigma'][new_index + max_index_offset]
    new_mean, new_sigma = add(old_mean, old_sigma, float(delay_mean), float(delay_sigma))
    new_arr['mean'][new_index + max_index_offset] = new_mean
    new_arr['sigma'][new_index + max_index_offset] = new_sigma

def add_arrival(arr, delay):
    new_arr = create_arrival()
    new_arr['pos_edge'] = arr['pos_edge']

    if delay['sense'] == "pos_unate":
        add_normal(arr, 0, delay, new_arr, 0)
        add_normal(arr, 1, delay, new_arr, 1)
    elif delay['sense'] == "neg_unate":
        add_normal(arr, 0, delay, new_arr, 1)
        add_normal(arr, 1, delay, new_arr, 0)
    elif delay['sense'] == "rising_edge":
        add_normal(arr, 0, delay, new_arr, 0)
        add_normal(arr, 0, delay, new_arr, 1)
    elif delay['sense'] == "falling_edge":
        add_normal(arr, 1, delay, new_arr, 0)
        add_normal(arr, 1, delay, new_arr, 1)
    else:
        print(f"ERROR: Unsupported sense {delay['sense']}")
    
    return new_arr

def merge_arrival(old_arr, new_arr):
    if old_arr == new_arr:
        return
    
    max_arrival(old_arr, new_arr, 0)
    max_arrival(old_arr, new_arr, 1)
    min_arrival(old_arr, new_arr, 2)
    min_arrival(old_arr, new_arr, 3)

def max_arrival(old_arr, new_arr, index):
    old_mean = old_arr['mean'][index]
    old_sigma = old_arr['sigma'][index]
    new_mean = new_arr['mean'][index]
    new_sigma = new_arr['sigma'][index]

    if new_mean == 1e99:
        new_arr['mean'][index] = 1e99
        new_arr['sigma'][index] = 0
    else:
        mean, sigma = max_normal(old_mean, old_sigma, new_mean, new_sigma)
        new_arr['mean'][index] = mean
        new_arr['sigma'][index] = sigma

def min_arrival(old_arr, new_arr, index):
    old_mean = old_arr['mean'][index]
    old_sigma = old_arr['sigma'][index]
    new_mean = new_arr['mean'][index]
    new_sigma = new_arr['sigma'][index]

    if new_mean == -1e99:
        new_arr['mean'][index] = -1e99
        new_arr['sigma'][index] = 0
    else:
        mean, sigma = min_normal(old_mean, old_sigma, new_mean, new_sigma)
        new_arr['mean'][index] = mean
        new_arr['sigma'][index] = sigma

def set_arrival_times(pin, arrs):
    for arr in arrs:
        add_arrival_time(pin, arr)


arrival_time = {}
def add_arrival_time(pin, arr):
    if pin not in arrival_time:
        arrival_time[pin] = []
    arrival_time[pin].append(arr)

def get_arrival_time(pin):
    return arrival_time.get(pin, [])


fanin={}
fanout={}

def add_fanin(from_node, to_node):
    if to_node not in fanin:
        fanin[to_node] = []
    fanin[to_node].append(from_node)

def add_fanout(from_node, to_node):
    if from_node not in fanout:
        fanout[from_node] = []
    fanout[from_node].append(to_node)

def get_fanin(pin):
    return fanin.get(pin, [])

def get_fanout(pin):
    return fanout.get(pin, [])

def get_all_vertices():
    all_pins = set()
    for from_node in fanout.keys():
        all_pins.add(from_node)
        to_pins = fanout[from_node]
        all_pins.update(to_pins)

    return list(all_pins)


def add_connection_data(line):
    from_node, to_node, *rest = map(str.strip, line.split(','))
    add_delay_arc(line)
    add_fanin(from_node, to_node)
    add_fanout(from_node, to_node)


check={}
def add_check_data(line):
    from_node, to_node, sense, rise_check, fall_check = map(str.strip, line.split(','))
    data = {'from': from_node, 'sense': sense, 'value': [rise_check, fall_check]}
    if to_node not in check:
        check[to_node]=[]
    check[to_node].append(data)

def get_check_data(to_node):
    return check.get(to_node,[])

# def read_graph_file(fname):
#     with open(fname, 'r') as file:
#         for line in file:
#             if file.tell() == 1 and line.startswith("#from vertex, to vertex"):
#                 continue
#             line = line.rstrip('\n')
#             add_connection_data(line)

def read_graph_file(fname):
    with open(fname, 'r') as file:
        first_line = next(file)  # Read the first line
        if not first_line.startswith("#from vertex, to vertex"):
            add_connection_data(first_line.strip())  # Process the first line

        for line in file:
            line = line.rstrip('\n')
            add_connection_data(line)

def init_fanin_counter():
    fanin_counter = {}
    all_pins = get_all_vertices()
    
    for pin in all_pins:
        fanin = get_fanin(pin)
        fanin_data = {'visited': 0, 'count': len(fanin)}
        fanin_counter[pin] = fanin_data

    return fanin_counter

#########################################################

def get_update_wavefront(fanin_counter):
    wavefront = []
    for pin, values in fanin_counter.items():
        if values.get('visited', 0) == 0 and values.get('count', 0) == 0:
            wavefront.append(pin)
            values['visited'] = 1
    return wavefront

def update_pin(pin):
    if debug>0:
        print(f"DEBUG: Calculating AT for {pin}")
    fanins = get_fanin(pin)

    arrivals = []
    for fanin in fanins:
        if debug > 1:
            print(f"DEBUG:   Calculating AT from {fanin}") 
        arcs = get_delay_arc(fanin, pin)
        fanin_arrs = get_arrival_time(fanin)

        for fanin_arr in fanin_arrs:
            if debug > 2:
                print(f"DEBUG:     Fanin AT  {'pos' if fanin_arr['pos_edge'] else 'neg'}_edge "
                      f"{{{fanin_arr['mean'][0]}:{fanin_arr['sigma'][0]} "
                      f"{fanin_arr['mean'][1]}:{fanin_arr['sigma'][1]} "
                      f"{fanin_arr['mean'][2]}:{fanin_arr['sigma'][2]} "
                      f"{fanin_arr['mean'][3]}:{fanin_arr['sigma'][3]}}} ")

            for arc in arcs:
                if debug > 3:
                    print(f"DEBUG:       Arc  {arc['sense']} "
                          f"{{{arc['mean'][0]}:{arc['sigma'][0]} "
                          f"{arc['mean'][1]}:{arc['sigma'][1]} "
                          f"{arc['mean'][2]}:{arc['sigma'][2]} "
                          f"{arc['mean'][3]}:{arc['sigma'][3]}}} ")

                new_arr = add_arrival(fanin_arr, arc)

                if debug > 3:
                    print(f"DEBUG:       New AT  {'pos' if new_arr['pos_edge'] else 'neg'}_edge "
                          f"{{{new_arr['mean'][0]}:{new_arr['sigma'][0]} "
                          f"{new_arr['mean'][1]}:{new_arr['sigma'][1]} "
                          f"{new_arr['mean'][2]}:{new_arr['sigma'][2]} "
                          f"{new_arr['mean'][3]}:{new_arr['sigma'][3]}}} ")

                if is_valid_arrival(new_arr):
                    arr_to_merge = find_match_arrival(new_arr, arrivals)
                    merge_arrival(new_arr, arr_to_merge)

                    if debug > 3:
                        print(f"DEBUG:       Final AT  {'pos' if arr_to_merge['pos_edge'] else 'neg'}_edge "
                              f"{{{arr_to_merge['mean'][0]}:{arr_to_merge['sigma'][0]} "
                              f"{arr_to_merge['mean'][1]}:{arr_to_merge['sigma'][1]} "
                              f"{arr_to_merge['mean'][2]}:{arr_to_merge['sigma'][2]} "
                              f"{arr_to_merge['mean'][3]}:{arr_to_merge['sigma'][3]}}} ")

            set_arrival_times(pin, arrivals)


def update_wavefront(wavefront):
    for pin in wavefront:
        update_pin(pin)  # Assuming update_pin is a function defined similarly to the Perl subroutine


def update_fanin_counter(wavefront, fanin_counter):
    for pin in wavefront:
        fanouts = get_fanout(pin)  # Assuming get_fanout is a function returning fanouts for a given pin
        for fanout in fanouts:
            if fanout in fanin_counter:
                fanin_counter[fanout]['count'] -= 1
                if fanin_counter[fanout]['count'] < 0:
                    print(f"ERROR: Negative fanin count in {fanout}")
                    exit()  # Exiting the program in case of an error


def calc_AT():
    fanin_counter = init_fanin_counter()  # Assuming init_fanin_counter initializes the fanin counter
    print(f"Total pins to calculate: {len(fanin_counter)}")

    wavefront = get_update_wavefront(fanin_counter)  # Assuming get_update_wavefront retrieves the wavefront
    init_root_arrival(wavefront)  # Assuming init_root_arrival initializes root arrival for wavefront

    iter_count = 0
    while len(wavefront) > 0:
        print(f"Iteration {iter_count}: {len(wavefront)} pins to calculate")

        update_wavefront(wavefront)  # Assuming update_wavefront updates the wavefront
        update_fanin_counter(wavefront, fanin_counter)  # Assuming update_fanin_counter updates fanin counters

        wavefront = get_update_wavefront(fanin_counter)  # Updating wavefront based on updated fanin counters
        iter_count += 1


def read_check_arc_data(fname):
    with open(fname, "r") as file:
        first_line = True
        for line in file:
            if first_line and "#from vertex, to vertex" in line:
                continue
            line = line.rstrip('\n')  # Removes newline character
            add_check_data(line)  # Assuming add_check_data processes the line
            first_line = False


def shifted_arrival(arr, index, N):
    if index < 2:
        return arr['mean'][index] + N * arr['sigma'][index]
    else:
        return arr['mean'][index] - N * arr['sigma'][index]


def calc_slack(data_arr, check_arr, check_arc):
    rise_slack = 1e99
    fall_slack = 1e99

    if check_arc['sense'] == "setup_rising":
        if check_arr['mean'][2] != 1e99:
            clock_shift = clock
            if check_arr['pos_edge'] != data_arr['pos_edge']:
                clock_shift = clock_shift / 2

            if data_arr['mean'][0] != -1e99:
                rise_slack = clock_shift + shifted_arrival(check_arr, 2, clock) - float(check_arc['value'][
                    0]) - float(shifted_arrival(data_arr, 0, clock))
                if debug > 0:
                    print(
                        f"DEBUG: {clock_shift} + {shifted_arrival(check_arr, 2, clock)} - {check_arc['value'][0]} - {shifted_arrival(data_arr, 0, clock)} = {rise_slack}")

            if data_arr['mean'][1] != -1e99:
                fall_slack = clock_shift + shifted_arrival(check_arr, 2, clock) - float(check_arc['value'][
                    1]) - float(shifted_arrival(data_arr, 1, clock))
                if debug > 0:
                    print(
                        f"DEBUG: {clock_shift} + {shifted_arrival(check_arr, 2, clock)} - {check_arc['value'][1]} - {shifted_arrival(data_arr, 1, clock)} = {fall_slack}")

    elif check_arc['sense'] == "setup_falling":
        if check_arr['mean'][3] != 1e99:
            clock_shift = clock
            if check_arr['pos_edge'] != data_arr['pos_edge']:
                clock_shift = clock_shift / 2

            if data_arr['mean'][0] != -1e99:
                rise_slack = clock_shift + shifted_arrival(check_arr, 3, clock) - check_arc['value'][
                    0] - shifted_arrival(data_arr, 0, clock)
                if debug > 0:
                    print(
                        f"DEBUG: {clock_shift} + {shifted_arrival(check_arr, 3, clock)} - {check_arc['value'][0]} - {shifted_arrival(data_arr, 0, clock)} = {rise_slack}")

            if data_arr['mean'][1] != -1e99:
                fall_slack = clock_shift + shifted_arrival(check_arr, 3, clock) - check_arc['value'][
                    1] - shifted_arrival(data_arr, 1, clock)
                if debug > 0:
                    print(
                        f"DEBUG: {clock_shift} + {shifted_arrival(check_arr, 3, clock)} - {check_arc['value'][1]} - {shifted_arrival(data_arr, 1, clock)} = {fall_slack}")

    else:
        print(f"ERROR: Unsupported sense {check_arc['sense']}")

    return rise_slack, fall_slack


def calc_end_slack(pin):
    rise_slack = 1e99
    fall_slack = 1e99

    check_arcs = get_check_data(pin)  # Assuming get_check_data retrieves check arcs for a pin
    data_arrs = get_arrival_time(pin)  # Assuming get_arrival_time retrieves arrival times for a pin

    for data_arr in data_arrs:
        for check_arc in check_arcs:
            check_pin = check_arc['from']
            check_arrs = get_arrival_time(
                check_pin)  # Assuming get_arrival_time retrieves arrival times for a check pin

            for check_arr in check_arrs:
                r, f = calc_slack(data_arr, check_arr, check_arc)
                rise_slack = min(rise_slack, r)
                fall_slack = min(fall_slack, f)

    return rise_slack, fall_slack


def calc_endpoint_slack(fname, out_fname):
    with open(fname, "r") as file_in, open(out_fname, "w") as file_out:
        for line in file_in:
            endpoint = line.strip()
            rise_slack, fall_slack = calc_end_slack(endpoint)
            file_out.write(f"{endpoint}, {rise_slack}, {fall_slack}\n")


def max_normal(mean_a, sigma_a, mean_b, sigma_b):
    shifted_a = mean_a + N * sigma_a
    shifted_b = mean_b + N * sigma_b
    if shifted_a > shifted_b:
        return mean_a, sigma_a
    else:
        return mean_b, sigma_b


def min_normal(mean_a, sigma_a, mean_b, sigma_b):
    shifted_a = mean_a + N * sigma_a
    shifted_b = mean_b + N * sigma_b
    if shifted_a < shifted_b:
        return mean_a, sigma_a
    else:
        return mean_b, sigma_b

def main():

    print("Read graph file ",graphCSVFile)
    read_graph_file(graphCSVFile)
    print("Calculating arrival time")
    calc_AT()
    print("Read check data file", checkFile)
    read_check_arc_data(checkFile)
    print("Calculate slack")
    calc_endpoint_slack(endpointsFile,outputSlackFile)
    print("Finished")


if __name__ == "__main__":
    main()