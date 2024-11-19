import yaml
import pickle

def load_yaml(filename):
    """
    Load a YAML file and return its contents as a Python dictionary,
    preserving the data types of values (int, float, list, str, etc.).
    
    Parameters
    ----------
    filename : str
        The path to the YAML file to load.
        
    Returns
    -------
    data : dict
        The contents of the YAML file as a dictionary with the correct data types.
    """
    try:
        with open(filename, 'r') as file:
            data = yaml.safe_load(file)  # Use safe_load to avoid loading any unsafe YAML tags
        return data
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file: {e}")

def save_pickle(data, filename):
    """
    Save data to a file using pickle.

    Parameters
    ----------
    data : obj 
        The data to be saved (can be any pickle-serializable object).
    filename : str
        The path to the file where data will be saved.

    Returns
    -------
    None
    """
    try:
        with open(filename, 'wb') as file:
            pickle.dump(data, file)
        print(f"Data successfully saved to '{filename}'")
    except Exception as e:
        print(f"Error saving data to pickle file: {e}")

def load_pickle(filename):
    """
    Load data from a pickle file.

    Parameters
    ----------
    filename : str
        The path to the pickle file to load.

    Returns
    -------
    data : obj
        The data loaded from the pickle file.
    """
    try:
        with open(filename, 'rb') as file:
            data = pickle.load(file)
        print(f"Data successfully loaded from '{filename}'")
        return data
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except pickle.PickleError as e:
        print(f"Error loading data from pickle file: {e}")
