import os


# -------------------------------
# STEP 1: Read full P-section data
# -------------------------------
def read_parameter_data(filename):
    param_data = []

    with open(filename, "r") as f:
        for line in f:
            section = line[72:73]

            if section == "P":
                data = line[:64]  # only IGES data field
                param_data.append(data)

    return param_data


# -------------------------------
# STEP 2: Reconstruct full records
# -------------------------------
def build_records(param_data):
    records = []
    current = ""

    for chunk in param_data:
        current += chunk.strip()

        if ";" in chunk:
            records.append(current)
            current = ""

    return records


# -------------------------------
# STEP 3: Parse type 110 entities
# -------------------------------
def parse_type110(records):
    coordinates = []

    for rec in records:
        rec = rec.replace(";", "")
        parts = rec.split(",")

        # Skip invalid or short records
        if len(parts) < 7:
            continue

        try:
            entity_type = int(parts[0])
        except:
            continue

        if entity_type == 110:
            try:
                nums = [float(x) for x in parts[1:7]]
                coordinates.append(nums)
            except ValueError:
                # Skip malformed entries safely
                continue


    return coordinates


# -------------------------------
# MAIN IMPORT FUNCTION
# -------------------------------
def import_SW(filename):
    param_data = read_parameter_data(filename)
    records = build_records(param_data)
    coords = parse_type110(records)

    return coords