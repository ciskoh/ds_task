""" 
    This contains the classes 
    that compose the Enzyme based model 
    presented by Allison et al. 2010
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List
import pathlib as p
from collections import Counter
import json

@dataclass
class AbstractProcess(ABC):

    outputs: Dict[str, float]
    
    @abstractmethod
    def calc(self, input):
        pass    

class EnzymeCatalysis(AbstractProcess):
    v_max: int =100000000
    km: int=500.0
    mic_to_soc: float = .5

    def calc(self, enz:float, soc:float,doc, eloss, mic, death=0):
        """ calculates DOC using Decomposition formula. WANRING: death is not calculated ans set equal to 0"""
        _assim = self.v_max*mic * (doc/(self.km+doc))
        _decomp = self.v_max*enz*(soc/(self.km+soc))
        ddoc= doc + ( death * (1-self.mic_to_soc) ) + _decomp + eloss - _assim
        self.outputs= {"doc":ddoc}

class Uptake(AbstractProcess):
    v_max_uptake: int= 100000000 
    km_uptake:float = 0.1
    
    def calc(self, mic, doc, cue, death=0, eprod=0):
        """calculates dMIC using microbial mass change formula 
        steps:
         - calculate assimilation
         - update cue calculation
         - calculater dmic
         - calculate co2 emissions 
        """
        _assim = self.v_max_uptake*mic * (doc/(self.km_uptake+doc))
        dmic = _assim *cue - death - eprod
        self.outputs = {"mic": dmic, "assim":_assim}

class CarbonUse(AbstractProcess):
    cue_slope:float = -0.016

    def calc(self, cue, temp, assim):
        _cue= self.cue_slope*temp + cue
        dco2 = assim*(1-_cue)
        self.outputs = {"co2": dco2}

class EnzymeProductionDecay(AbstractProcess):
    r_enz_prod: float = 0.000005
    r_enz_loss: float = 0.001
    def calc(self, mic, enz):
        _eprod = mic*self.r_enz_prod
        _eloss = enz * self.r_enz_loss
        
        d_enz = _eprod - _eloss
        self.outputs={"enz": d_enz, "eprod":_eprod, "eloss": _eloss}

class MicDeath(AbstractProcess):
    r_death: float= 0.0002 

    def calc(self,mic:float):
        death=self.r_death *mic
        self.outputs = {"death": death}


@dataclass
class ModelRunner:
    version: str = 0.1
    #main variables
    soc: float = None
    doc: float = None
    mic: float = None
    enz: float = None
    # extra intermediate variables, used for computation
    eloss:float = None
    death:float=None
    cue:float=None
    eprod: float=None
    assim:float=None
    # output variables
    co2: float = None
    # processes defined above
    encat: EnzymeCatalysis
    uptake: Uptake
    cbuse: CarbonUse
    enzpd: EnzymeProductionDecay
    micdeath: MicDeath

    history: List[Dict[str, float]]

    def compute_model(self, time, temp):
        
        print("computing model for dtime={dtime}, dtemp={dtemp}")
        if not self.eloss:
            self.eloss= self.enz*self.enzpd.r_enz_loss
        if not self.eprod:
            self.eprod=self.mic * self.enzpd.r_enz_prod
        if not self.death:
            self.death=self.micdeath.r_death*self.mic
        
        self.encat.calc(enz=self.enz, soc=self.soc, mic=self.mic, death=self.death)
        self.uptake.calc(temp, self.mic, self.doc, self.death,self.eprod)
        self.cbuse.calc(self.cue, temp, self.uptake.outputs["assim"] )
        self.enzpd.calc(self.mic, self.enz)
        self.micdeath.calc(self.mic)        
        print("computation finished")
        self.collect_results()
        self.save_history()
        print(f"current state of the system:{self.history[-1]=}")
              
    def collect_results(self):
        res = Counter(self.encat.outputs).update(self.uptake.outputs).update(self.cbuse.outputs)
        
        self.soc+=res.get("soc",0)
        self.doc+=res.get("doc",0)
        self.mic+=res.get("mic",0)
        self.enz+=res.get("enz",0)
        self.co2+=res.get("co2",0)

    def save_history(self):
        self.history.append(self.to_dict())
    
    def export_history(self, dst_path: str | p.Path):
        with open(dst_path, "w+") as f:
            json.dump(self.history, f)
            if p.Path(dst_path).is_file() and p.Path(dst_path).exists():
                print(f"computation history saved at {dst_path}")
            else:
                raise FileNotFoundError(f"could not find history file at {dst_path}. Something went wrong")

