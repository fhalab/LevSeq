import os

def check_if_folder_exists(folder_path):
    """ Check if a folder exists in the results folder exists. If not, create a new one.
    Input: folder_path
    Output: True if folder exists, False otherwise
    """
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        


