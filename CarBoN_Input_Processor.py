import numpy as np
import pandas as pd


class Grains:
    def __init__(self,r_min,r_max,Ratio,Hamaker,density,num):
        self.column_names = ['Input1','Input2','Output1','Output2',
                        'f_ijk','K_ij','Hamaker','Radius1','Radius2','Fo']
        self.r_min = r_min
        self.r_max = r_max
        self.R = Ratio    
        self.A = Hamaker
        self.p = density
        self.rows_list = []
        self.num = num
    def _create_bins(self): 
        self.numbins = int(np.ceil(1+np.log((self.r_max/self.r_min)**3)/np.log(self.R)))
        self.r,self.v,self.m=np.zeros(self.numbins+1),np.zeros(self.numbins+1),\
                              np.zeros(self.numbins+1)
        for n in np.arange(1,self.numbins+1):
            self.r[n] = self.r_min*self.R**((n-1)/3)
            self.v[n] = 4/3*np.pi*self.r_min**3*self.R**(n-1)
            self.m[n] = self.p*4/3*np.pi*(self.r_min*100)**3*self.R**(n-1)
    def _create_coeffs(self,r,v,m,i,j,R):
        k = 1.380649e-23 # Boltzmann Constant
        rij = (r[i]**3+r[j]**3)**(1/3)
        vij = v[i]+v[j]
        muij = (m[i]*m[j])/(m[i]+m[j])
        Kij = np.pi*(r[i]+r[j])**2*np.sqrt((8*k)/(np.pi*muij))
        bin_num = int(np.floor(1+3*np.log(rij/r[1])/np.log(R)))
        fraction = (r[bin_num+1]**3-rij**3)/(r[bin_num+1]**3-r[bin_num]**3)
        fraction2 = (r[bin_num+1]**3-rij**3)/(r[bin_num+1]**3-r[bin_num]**3)*(r[bin_num]/rij)
        return rij,vij,Kij,bin_num,fraction,fraction2
    def _create_reac_df(self):
        self._create_bins()
        for i in np.arange(1,self.numbins):
            for j in np.arange(1,i+1):
                rij,vij,Kij,bin_num,fraction,fraction2=self._create_coeffs(self.r,self.v,\
                                                                           self.m,i,j,self.R)
                reac_coeffs=[i+self.num,j+self.num,bin_num+self.num,bin_num+self.num+1,\
                             fraction2,Kij,self.A,self.r[i],self.r[j],6]
                reac_row={k:v for k,v in zip(self.column_names,reac_coeffs)}
                self.rows_list.append(reac_row)
    def output(self):
        self._create_reac_df()
        return pd.DataFrame(self.rows_list,columns=self.column_names)

#

class Kida:
    def __init__(self,reac_file,spec_file):
        self.reac_col_names = ['Input1','Input2','Output1','Output2','Output3',
                        'alpha','beta','gamma','F','g','Type','Re','Tlo', 
                        'Thi','Fo','N','V','R']
        self.reac_col_widths=[11,23,11,10,34,11,11,11,9,9,5,3,7,7,3,5,2,3]
        self.reac_dtypes={'Input1':str,'Input2':str,'Output1':str,'Output2':str,
                          'Output3':str}
        self.reac_file = reac_file
        self.spec_file = spec_file

    def read_reactions(self):
        if not hasattr(self, '_spec_dict'):
            raise Exception("Must read species first")

        self._reac_df = pd.read_fwf(self.reac_file, comment='#', header=None,
                          names=self.reac_col_names,widths=self.reac_col_widths,
                          converters=self.reac_dtypes)
        self._process_reactions_file()

    def read_species(self):
        self._spec_df = pd.read_fwf(self.spec_file, comment='#', header=None)
        #print(self._spec_df)
        self._process_species_file()

    def species_dictionary(self):
        return self._spec_dict

    def species_dataframe(self):
        return self._spec_df

    def reactions_dataframe(self):
        return self._reac_df

    def output(self):
        return  self._reac_df,self._spec_df, self._spec_dict

    def _process_species_file(self):
        col_list=list(self._spec_df)
        self._spec_df['atom_num'] = self._spec_df[col_list[2:24]].sum(axis=1)
        self._spec_df = self._spec_df.drop(col_list[2:24], axis=1)
        self._spec_df.rename(columns={0 : "species",
                           1 : "charge",
                           24: "species_num"}, inplace=True)

        self._spec_dict = pd.Series(self._spec_df.species_num.values,
                               index=self._spec_df.species).to_dict()
        add_photons={'Photon':0,'Pho':0}
        self._spec_dict.update(add_photons)
        self.num_species=max(self._spec_df.species_num.values) 

    def _process_reactions_file(self):
        for column in self._reac_df.columns[0:5]: 
            self._reac_df[column] = self._reac_df[column].replace(self._spec_dict)
            self._reac_df[column] = self._reac_df[column].astype('Int64')


reac = "data/kida_reac_C_O_only.dat"
spec = "data/kida_spec_C_O_only.dat"

K=Kida(reac,spec)

K.read_species()
K.read_reactions()
reactions, species, dictionary = K.output()

print(reactions)
print(species)
print(dictionary)

g=Grains(1e-11,1e-6,2,2e-20,2.3,K.num_species)



grain_reacs = g.output()

print(grain_reacs)


