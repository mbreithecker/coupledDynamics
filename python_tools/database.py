"""
Helper script to load and filter simulation data created by the cuda/c++ program
"""
import json


class Database:
    DATA = list()

    def __init__(self, data=None):
        if data is None:
            data = list()

        self.DATA = data

    def load_file(self, file):
        """
        Loads all entries of a simulation file and adds it to the database.
        All parameters named in the file are preserved.
        :param file: Name of a measurement-file generated by the cpp program
        :return: none
        """
        with open(file) as f:
            self.DATA += json.loads(
                ("[" + f.read() + "]").replace("\n", " ").replace(", ]", "]").replace(",]", "]").replace(",}", "}"))

        return self

    def get_data_grouped_by(self, key_list):
        """
        Will return a list of a lists with measure_entries,
        On the first layer, every measure_entry has the same value for every key specified in key-list
        :param key_list: parameters that should be the same for a entry group
        :return: list of lists, all entries in the nested lists are equals for the keys in key_list
        """
        groups = dict()

        for entry in self.DATA:
            key_string = ""
            for key in key_list:
                key_string += str(entry.get(key)) + "-%-"

            if key_string not in groups:
                groups[key_string] = list()

            groups[key_string].append(entry)

        return list(groups.values())

    def filter_by_parameters(self, param_value_dict):
        filtered_entries = list()
        for entry in self.DATA:
            include_flag = True
            for param, value in param_value_dict.items():
                if entry.get(param) != value:
                    include_flag = False
                    break

            if include_flag:
                filtered_entries.append(entry)

        return Database(filtered_entries)

    def __print_statistics(self, entries, excluded_keys, depth, ignore_keys_list):
        key_list = list()
        remaining_groups_list = list()
        value_tree = dict()
        for entry in entries:
            for sim_key in entry.keys():
                sim_value = entry.get(sim_key)
                if type(sim_value) is list:
                    continue

                if sim_key in ignore_keys_list:
                    continue

                if sim_key not in value_tree:
                    value_tree[sim_key] = dict()

                if sim_value not in value_tree[sim_key]:
                    value_tree[sim_key][sim_value] = 0

                value_tree[sim_key][sim_value] += 1

        # count same
        lowest_key_name = ""
        lowest_key_count = 10e10
        for key in value_tree.keys():
            if len(value_tree.get(key).values()) == 1:
                key_list.append(key)
            else:
                entry_len = len(value_tree.get(key).keys())
                if entry_len < lowest_key_count:
                    lowest_key_count = entry_len
                    lowest_key_name = key

        if lowest_key_name != "":
            for val in value_tree.get(lowest_key_name).keys():
                group = list()
                for entry in entries:
                    if val == entry.get(lowest_key_name):
                        group.append(entry)

                remaining_groups_list.append(group)

        print(" -" * depth, [[p, entries[0].get(p)] for p in key_list if p not in excluded_keys], " amount: ",
              len(entries))
        for entry_group in remaining_groups_list:
            self.__print_statistics(entry_group, excluded_keys + key_list, depth + 1, ignore_keys_list)

    def print_statistics(self, ignore_keys_list=None):
        if ignore_keys_list is None:
            ignore_keys_list = list()
        self.__print_statistics(self.DATA, [], 0, ignore_keys_list)

    def merge_measurements_with_same_parameters(self):
        key_set = set()
        merge_keys = set()
        for entry in self.DATA:
            for key, val in entry.items():
                if type(val) != list and key != "seed":
                    key_set.add(key)
                elif type(val) == list:
                    merge_keys.add(key)

        groups = self.get_data_grouped_by(list(key_set))

        merged_list = list()
        for group in groups:
            # Merge
            merged_entry = group[0]
            for key in merge_keys:
                for i in range(1, len(group)):
                    merged_entry[key] += group[i][key]

            merged_entry["A"] = len(merged_entry.get("sigma"))
            merged_list.append(merged_entry)

        self.DATA = merged_list
