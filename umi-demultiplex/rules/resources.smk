def get_mem_20_10(wildcards, attempt):
    return str(10 + attempt * 10) + "G"

def get_mem_3_1(wildcards, attempt):
    return str(2 + attempt) + "G"

def get_mem_30_10(wildcards, attempt):
    return str(20 + attempt * 10) + "G"

def get_mem_40_10(wildcards, attempt):
    return str(30 + attempt * 10) + "G"

def get_mem_70_10(wildcards, attempt):
    return str(60 + attempt * 10) + "G"

def get_mem_100_10(wildcards, attempt):
    return str(90 + attempt * 10) + "G"

def get_mem_100_50(wildcards, attempt):
    return str(90 + attempt * 50) + "G"

def get_time_2_1(wildcards, attempt):
    return "0" + str(1 + attempt) + ":00:00"

def get_time_3_1(wildcards, attempt):
    return "0" + str(2 + attempt) + ":00:00"

def get_time_1_1(wildcards, attempt):
    return "0" + str(0 + attempt) + ":00:00"
