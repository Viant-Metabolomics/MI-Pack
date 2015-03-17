#!/usr/local/bin/python

import sqlite3

class SQLiteManager(object):
    def __init__(self, dbname):
        self.connection = sqlite3.connect('%s' % dbname)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()

    def create_table(self, table_type, ion_mode, name_peak_list = "", sql_str = ""):

        self.name_peak_list = name_peak_list

        if table_type == "SPS":
            self.cursor.execute("DROP TABLE IF EXISTS SPS_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE SPS_%s_%s (
                KEGG_COMPOUND char(6) DEFAULT NULL,
                ion char(20) DEFAULT NULL,
                formula text DEFAULT NULL,
                name text DEFAULT NULL,
                mz decimal(12,7) DEFAULT NULL,
                intensity decimal(18,8) DEFAULT NULL,
                mass decimal(12,7) DEFAULT NULL,
                ppm_error decimal(12,2) DEFAULT NULL,
                exact_mass decimal(12,7) DEFAULT NULL,
                primary key  (KEGG_COMPOUND, ion, mz, intensity)
                );""" % (ion_mode, self.name_peak_list))

        if table_type == "TM":
            self.cursor.execute("DROP TABLE IF EXISTS transformations_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE transformations_%s_%s (
                mz_a decimal(12,7) DEFAULT NULL,
                mz_b decimal(12,7) DEFAULT NULL,
                ion_a char(20) DEFAULT NULL,
                ion_b char(20) DEFAULT NULL,
                intensity_a decimal(18,8) DEFAULT NULL,
                intensity_b decimal(18,8) DEFAULT NULL,        
                transformation char(50) DEFAULT NULL,
                mz_pair_diff decimal(12,7) DEFAULT NULL,
                mz_pair_diff_neutral decimal(12,7) DEFAULT NULL,
                mass_pair_diff decimal(12,7) DEFAULT NULL,
                ppm_error decimal(12,2) DEFAULT NULL,
                type int(2) DEFAULT NULL
                );""" % (ion_mode, self.name_peak_list))


        if table_type == "compounds":

            self.cursor.execute("""DROP TABLE IF EXISTS compounds_%s;""" % (ion_mode))
            self.cursor.execute("""DROP TABLE IF EXISTS map_compound;""")
            self.cursor.execute("""DROP TABLE IF EXISTS map_names;""")
            self.cursor.execute("""DROP TABLE IF EXISTS map_formula;""")

            self.cursor.execute("""
                CREATE TABLE compounds_%s (
                KEGG_COMPOUND char(6) not null,
                name text,
                formula text,
                exact_mass decimal(12,7) default null,
                mass decimal(12,7) default null,
                ion char(20),
                PRIMARY KEY (KEGG_COMPOUND, ion)
                );""" % (ion_mode))

            self.cursor.execute("""
                CREATE TABLE map_names (
                mapid char(8) not null,
                mapname tinytext,
                PRIMARY KEY (mapid)
                );""")

            self.cursor.execute("""
                CREATE TABLE map_compound (
                KEGG_COMPOUND char(6) not null,
                mapid char(8) not null,
                PRIMARY KEY (mapid, KEGG_COMPOUND)
                );""")

            self.cursor.execute("""
                CREATE TABLE map_formula (
                mapid char(8) not null,
                formula text not null,
                PRIMARY KEY (mapid, formula)
                );""")

        if table_type == "reactant_pairs":
                
            self.cursor.execute("drop TABLE if exists reactions;")
            self.cursor.execute("drop TABLE if exists reactant_pairs_%s;" % (ion_mode))
                
            self.cursor.execute("""
                CREATE TABLE reactions (
                id mediumint(9) not null,
                reactionid char(6) not null,
                type char(1) default null,
                PRIMARY KEY (id, reactionid, type)
                );""")

            self.cursor.execute("""
                CREATE TABLE reactant_pairs_%s (
                id int(11) default null,
                KEGG_COMPOUND_a char(6) not null,
                KEGG_COMPOUND_b char(6) not null,
                ion_a char(20) default null,
                ion_b char(20) default null,
                type int(1) default null,
                mass_pair_diff decimal(12,7) default null,
                mz_pair_diff decimal(12,7) default null,
                transformation char(50) default null,
                PRIMARY KEY  (KEGG_COMPOUND_a, KEGG_COMPOUND_b, ion_a, ion_b, type)
                );""" % (ion_mode))

        if table_type == "organism":
            self.cursor.execute("""
                CREATE TABLE organism_compound (
                organismid char(3) not null,
                KEGG_COMPOUND char(6) not null,
                PRIMARY KEY (organismid, KEGG_COMPOUND)
                );""")

        if table_type == "temp_reactant_pairs":
            self.cursor.execute("DROP TABLE IF EXISTS temp_reactant_pairs;")
            self.cursor.execute("""
                    CREATE TEMPORARY TABLE temp_reactant_pairs (
                      id1 int(11) not null,
                      KEGG_COMPOUND_a char(6) not null,
                      id2 int(11) not null,
                      KEGG_COMPOUND_b char(6) not null
                    );
                    """)
            self.cursor.execute("""INSERT INTO temp_reactant_pairs SELECT DISTINCT RPS.id AS id1, RPS.KEGG_COMPOUND_a, cs.id AS id2, cs.KEGG_COMPOUND_b
                FROM reactant_pairs_%s as RPS inner join 
                reactant_pairs_%s AS cs on RPS.KEGG_COMPOUND_b = cs.KEGG_COMPOUND_a and RPS.KEGG_COMPOUND_a != cs.KEGG_COMPOUND_b;""" % (ion_mode, ion_mode))
            self.save()

        if table_type == "unique_transformations":
            self.cursor.execute("DROP TABLE IF EXISTS unique_transformations;")

            self.cursor.execute("""CREATE TABLE unique_transformations AS SELECT DISTINCT
                transformation, mass_pair_diff, type
                FROM rpairs ORDER BY mass_pair_diff;""")
            
            self.cursor.execute("""CREATE UNIQUE INDEX idx_unique_transformations ON unique_transformations
                (transformation, mass_pair_diff, type);""")

        if table_type == "SUBSET_unqiue_transformations":
                       
            self.cursor.execute("DROP TABLE IF EXISTS SUBSET_unique_transformations;")
            
            self.cursor.execute("""CREATE TEMPORARY TABLE SUBSET_unique_transformations AS SELECT DISTINCT
                transformation, mass_pair_diff, type
                FROM SUBSET_rpairs ORDER BY mass_pair_diff;""")

            self.cursor.execute("""DROP INDEX IF EXISTS PK_SUBSET_unique_transformations""")  
            self.cursor.execute("""CREATE UNIQUE INDEX PK_SUBSET_unique_transformations ON SUBSET_unique_transformations
                (transformation, mass_pair_diff, type);""")

            
        if table_type == "reactant_pairs_peaklist":
            self.cursor.execute("DROP TABLE IF EXISTS reactant_pairs_%s_%s_temp;" % (ion_mode, self.name_peak_list))
            #TEMPORARY 
            self.cursor.execute("""CREATE TABLE reactant_pairs_%s_%s_temp (
                KEGG_COMPOUND_a char(6) DEFAULT NULL,
                ion_a char(20) DEFAULT NULL,
                KEGG_COMPOUND_b char(6) DEFAULT NULL,
                ion_b char(20) DEFAULT NULL,
                transformation char(50) DEFAULT NULL,
                ppm_error decimal(12,2) DEFAULT NULL,        
                mz_a decimal(12,7) DEFAULT NULL,
                mz_b decimal(12,7) DEFAULT NULL,
                intensity_a decimal(18,8) DEFAULT NULL,
                intensity_b decimal(18,8) DEFAULT NULL, 
                mz_pair_diff decimal(12,7) DEFAULT NULL,
                type int(1) DEFAULT NULL
                );""" % (ion_mode, self.name_peak_list))

            self.cursor.execute("DROP TABLE IF EXISTS reactant_pairs_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""CREATE TABLE reactant_pairs_%s_%s (
                KEGG_COMPOUND_a char(6) DEFAULT NULL,
                ion_a char(20) DEFAULT NULL,
                KEGG_COMPOUND_b char(6) DEFAULT NULL,
                ion_b char(20) DEFAULT NULL,
                transformation char(50) DEFAULT NULL,
                ppm_error decimal(12,2) DEFAULT NULL,        
                mz_a decimal(12,7) DEFAULT NULL,
                mz_b decimal(12,7) DEFAULT NULL,
                intensity_a decimal(18,8) DEFAULT NULL,
                intensity_b decimal(18,8) DEFAULT NULL,
                mz_pair_diff decimal(12,7) DEFAULT NULL,
                type int(1) DEFAULT NULL
                );""" % (ion_mode, self.name_peak_list))

            self.cursor.execute("DROP TABLE IF EXISTS TM_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""CREATE TABLE TM_%s_%s (
                KEGG_COMPOUND char(6) NOT NULL,
                ion char(20) DEFAULT NULL,
                formula text DEFAULT NULL,
                name text DEFAULT NULL,
                mz decimal(12,7) DEFAULT NULL,
                intensity decimal(18,8) DEFAULT NULL,
                mass decimal(12,7) DEFAULT NULL,
                ppm_error decimal(12,2) DEFAULT NULL,
                exact_mass decimal(12,7) DEFAULT NULL,
                type int(1) DEFAULT NULL,
                primary key  (KEGG_COMPOUND, ion, mz, type)
                );""" % (ion_mode, self.name_peak_list))
            
        if table_type == "maps_SPS":
            self.cursor.execute("DROP TABLE IF EXISTS map_cpd_SPS_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE map_cpd_SPS_%s_%s (
                mapid char(8) not null,
                KEGG_COMPOUND char(6) NOT NULL,
                primary key  (KEGG_COMPOUND, mapid)
                );""" % (ion_mode, self.name_peak_list))

            self.cursor.execute("DROP TABLE IF EXISTS map_f_SPS_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE map_f_SPS_%s_%s (
                mapid char(8) not null,
                formula text NOT NULL,
                primary key  (formula, mapid)
                );""" % (ion_mode, self.name_peak_list))

            self.cursor.execute("DROP TABLE IF EXISTS maps_coverage_SPS_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE maps_coverage_SPS_%s_%s (
                mapid char(8) not null,
                mapname text not null,
                url text NOT NULL,
                subtotal int(11) NOT NULL,
                total int(11) NOT NULL,
                coverage int(11) NOT NULL,
                identity char(1) NOT NULL,
                primary key (mapid, identity)
                );""" % (ion_mode, self.name_peak_list))
            
        if table_type == "maps_TM":

            self.cursor.execute("DROP TABLE IF EXISTS map_cpd_TM_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE map_cpd_TM_%s_%s (
                mapid char(8) not null,
                KEGG_COMPOUND char(6) NOT NULL,
                primary key  (KEGG_COMPOUND, mapid)
                );""" % (ion_mode, self.name_peak_list))

            self.cursor.execute("DROP TABLE IF EXISTS map_f_TM_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE map_f_TM_%s_%s (
                mapid char(8) not null,
                formula text NOT NULL,
                primary key  (formula, mapid)
                );""" % (ion_mode, self.name_peak_list))

            self.cursor.execute("DROP TABLE IF EXISTS maps_coverage_TM_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE maps_coverage_TM_%s_%s (
                mapid char(8) not null,
                mapname text not null,
                url text NOT NULL,
                subtotal int(11) NOT NULL,
                total int(11) NOT NULL,
                coverage int(11) NOT NULL,
                identity char(1) NOT NULL,
                primary key (mapid, identity)
                );""" % (ion_mode, self.name_peak_list))
                     
        if table_type == "formulae":
            self.cursor.execute("DROP TABLE IF EXISTS EFS_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""CREATE TABLE EFS_%s_%s (
                mz decimal(12,7) DEFAULT NULL,
                intensity decimal(16,8) DEFAULT NULL,
                mass decimal(12,7) DEFAULT NULL,
                exact_mass decimal(12,7) DEFAULT NULL,
                ppm_error decimal(12,7) DEFAULT NULL,
                formula text DEFAULT NULL,
                ion char(15) DEFAULT NULL,
                lewis int(1) DEFAULT NULL,
                senior int(1) DEFAULT NULL,
                HC int(1) DEFAULT NULL,
                NOPSC int(1) DEFAULT NULL,
                rules int(1) DEFAULT NULL,
                PRIMARY KEY (mz, intensity, formula, ion)
                );""" % (ion_mode, self.name_peak_list))

        if table_type == "peak_patterns":
            self.cursor.execute("DROP TABLE IF EXISTS PPS_%s_%s;" % (ion_mode, self.name_peak_list))
            self.cursor.execute("""
                CREATE TABLE PPS_%s_%s (
                id int(11) DEFAULT NULL,
                mz_pair_diff decimal(12,7) DEFAULT NULL,
                mz_a decimal(12,7) DEFAULT NULL,
                mz_b decimal(12,7) DEFAULT NULL,
                intensity_a decimal(18,8) DEFAULT NULL,
                intensity_b decimal(18,8) DEFAULT NULL,
                SNR_mean_a decimal(18,8) DEFAULT NULL,
                SNR_mean_b decimal(18,8) DEFAULT NULL,
                label_a char(15) DEFAULT NULL,
                label_b char(15) DEFAULT NULL,
                type text DEFAULT NULL,
                num_peaks int(2) DEFAULT NULL,
                n_atoms int DEFAULT NULL,
                mean_atoms decimal(12,7) DEFAULT NULL,
                sd_atoms decimal(12,7) DEFAULT NULL,
                rules int(1) DEFAULT NULL,
                PRIMARY KEY (id, mz_a, mz_b, label_a, label_b)
                );""" % (ion_mode, self.name_peak_list))

    def execute_complex(self, query_type, input_list, name_peak_list = "", ion_mode="", subset=""):

        if name_peak_list != "":
            self.name_peak_list = name_peak_list

        if query_type == "reactant_pairs_temp":
            
            self.cursor.execute("""
            INSERT INTO reactant_pairs_%s_%s_temp
            SELECT RP.KEGG_COMPOUND_a, diffs.ion_a, RP.KEGG_COMPOUND_b, diffs.ion_b,
            RP.transformation, ppm_error, mz_a, mz_b, intensity_a, intensity_b, diffs.mz_pair_diff,
            RP.type

            FROM %srpairs as RP INNER JOIN transformations_%s_%s as diffs WHERE
        
            RP.KEGG_COMPOUND_a = '%s' AND diffs.mz_a = %s
            AND diffs.ion_a = '%s'
            AND RP.transformation = diffs.transformation
            AND RP.type <= %s AND RP.type = diffs.type
            
            """ % tuple(input_list))
            # AND RP.ion_mode = '%s' AND RP.ion_mode = diffs.ion_mode;

            self.save()
            
        if query_type == "maps_TM":
            
            self.cursor.execute("""INSERT INTO map_cpd_TM_%s_%s (mapid, KEGG_COMPOUND)
                SELECT distinct map_compound.mapid, TM.KEGG_COMPOUND
                FROM TM_%s_%s as TM
                INNER JOIN map_compound WHERE
                TM.KEGG_COMPOUND = map_compound.KEGG_COMPOUND 
                %s
                ;""" % (ion_mode, self.name_peak_list, ion_mode, self.name_peak_list, input_list[0]))

            self.cursor.execute("""INSERT INTO map_f_TM_%s_%s (mapid, formula)
                SELECT distinct map_formula.mapid, TM.formula
                FROM TM_%s_%s as TM
                INNER JOIN map_formula WHERE
                TM.formula = map_formula.formula
                %s
                ;""" % (ion_mode, self.name_peak_list, ion_mode, self.name_peak_list, input_list[1]))

            self.save()

        if query_type == "maps_SPS":
            self.cursor.execute("""INSERT INTO map_cpd_SPS_%s_%s (mapid, KEGG_COMPOUND)
                SELECT distinct map_compound.mapid, SPS.KEGG_COMPOUND
                FROM SPS_%s_%s as SPS
                INNER JOIN map_compound WHERE
                SPS.KEGG_COMPOUND = map_compound.KEGG_COMPOUND
                %s
                ;""" % (ion_mode, self.name_peak_list, ion_mode, self.name_peak_list, input_list[0]))

            self.cursor.execute("""INSERT INTO map_f_SPS_%s_%s (mapid, formula)
                SELECT distinct map_formula.mapid, SPS.formula
                FROM SPS_%s_%s as SPS
                INNER JOIN map_formula WHERE
                SPS.formula = map_formula.formula
                %s
                ;""" % (ion_mode, self.name_peak_list, ion_mode, self.name_peak_list, input_list[1]))
            self.save()
            
    def index(self, index_name, table_name, columns):

        if len(columns) > 1:
            columns = ", ".join(columns)
        else:
            columns = columns[0]
        
        self.cursor.execute("""DROP INDEX IF EXISTS idx_%s""" % (index_name))  
        self.cursor.execute("""
        CREATE INDEX idx_%s
        on %s (%s);
        """ % (index_name, table_name, columns))

    def attach(self, ion_mode, db_name_to_attach, path_MIDB):
        self.cursor.execute("attach '%s' as '%s'" % (path_MIDB, db_name_to_attach))
       
    def drop(self, table_name):
        self.cursor.execute("DROP TABLE IF EXISTS %s;" % (table_name))
       
    def __len__(self):
        a = 0
        for i in self.__iter__():
            a += 1
        return a
       
    def execute(self, command):
        self.cursor.execute(command)
                  
    def max(self, table_name, column, where=""):
        self.cursor.execute("SELECT MAX(%s) FROM %s %s;" % (column, table_name, where))
        return self.cursor.fetchone()[0]

    def min(self, table_name, column, where=""):
        self.cursor.execute("SELECT MIN(%s) FROM %s %s;" % (column, table_name, where))
        return self.cursor.fetchone()[0]

    def count(self, table_name, column, where=""):
        self.cursor.execute("SELECT count(%s) FROM %s %s;" % (column, table_name, where))
        return self.cursor.fetchone()[0]
       
    def select(self, table_name, columns, where=""):
        if len(columns) > 1:
            columns = ", ".join(map(str,columns))
        else:
               columns = columns[0]       
        self.cursor.execute("SELECT %s FROM %s %s" % (columns, table_name, where))
        return self.cursor.fetchall()
                   
    def insert(self, table_name, columns_values):
        self.cursor.execute("INSERT INTO %s (%s) VALUES (%s)" % (table_name, ", ".join(map(str,columns_values.keys())), str(columns_values.values())[1 : -1]))
        
    def delete(self, table_name, where=""):
        self.cursor.execute("DELETE FROM %s %s" % (table_name, where))

    def create_select(self, table_name_insert, columns_select, table_name_select, where="", temporary=0):

        if temporary == 1:
            temporary = "TEMPORARY"
        else:
            temporary = ""
        
        if len(columns_select) > 1:
            columns_select = ", ".join(map(str,columns_select))
        else:
            columns_select = columns_select[0]

        self.cursor.execute("""CREATE %s TABLE %s as SELECT %s FROM %s %s;
        """ % (temporary, table_name_insert, columns_select, table_name_select, where))


    def insert_select(self, table_name_insert, columns_insert, columns_select, table_name_select, where=""):

        if len(columns_insert) > 1:
            columns_insert = ", ".join(map(str,columns_insert))
        else:
            columns_insert = columns_insert[0]

        if len(columns_select) > 1:
            columns_select = ", ".join(map(str,columns_select))
        else:
            columns_select = columns[0]

        self.cursor.execute("""INSERT INTO %s (%s)
        SELECT distinct %s FROM %s %s;
        """ % (table_name_insert, columns_insert, columns_select, table_name_select, where))
           
    def fetchall(self, values_only=""):
        if values_only != "":
            trueValues = []
            for value in self.cursor.fetchall():
                if len(value) > 1:
                    trueValues.append(value)
                else:
                    trueValues.append(value[0])
            return trueValues
        else:
            return self.cursor.fetchall()
       
    def fetchone(self):
        value = self.cursor.fetchone()
        if value != None:
            return value[0]
        else:
            return None

    def save(self):
        self.connection.commit()
        
    def close(self):
        self.connection.close()

