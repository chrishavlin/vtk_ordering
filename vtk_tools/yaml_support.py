import yaml
import os
class folded_str(str): pass
class literal_str(str): pass

def folded_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='>')
def literal_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(folded_str, folded_str_representer)
yaml.add_representer(literal_str, literal_str_representer)


class yt_helper(object):
    def __init__(self):
        try: 
            import yt 
        except ImportError:
            raise ImportError("yt_helper requires yt! Failed to import yt.")
            
        self.mesh_yaml = self._get_yt_mesh_yaml_path()
        
        with open(self.mesh_yaml,'r') as yml:
            self.mesh_types = yaml.full_load(yml)
        
        self.available_mesh_types=list(self.mesh_types.keys())
        
    def _get_yt_mesh_yaml_path(self):
        try: 
            import yt 
        except ImportError:
            raise ImportError("get_yt_mesh_yaml_path() requires yt! Failed to import yt.")

        yt_path =yt.__path__[0]
        mesh_yaml = os.path.join(yt_path,'utilities','mesh_types.yaml')
        if os.path.isfile(mesh_yaml):
            return mesh_yaml 
        else:
            raise FileNotFoundError(f"Could not find {mesh_yaml}")
            
    
