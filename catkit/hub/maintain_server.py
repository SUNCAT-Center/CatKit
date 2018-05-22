from catkit.hub.postgresql import CathubPostgreSQL


class MaintainPostgres(CathubPostgreSQL):
    def fill_reaction_system(self):
        con = self.connection or self._connect()
        cur = con.cursor()
        cur.execute(
            'SELECT distinct id from reaction where id not in (SELECT distinct id from reaction_system);')
        result = cur.fetchall()
        for id in result:
            id = id[0]
            cur.execute(
                "INSERT INTO reaction_system(name, id, ase_id) VALUES ('N/A', {}, '214d69b1b872fbfd5bb017e05153eaaa');".format(id))
        con.commit()
        con.close()

    def delete_lost_systems(self):
        con = self.connection or self._connect()
        cur = con.cursor()
        cur.execute(
            'SELECT distinct id from systems where unique_id not in (SELECT distinct ase_id from reaction_system);')
        result = cur.fetchall()

        for table in ['keys, text_key_values, number_key_values, species']:
            cur.execute(
                "delete from {} where id in (SELECT distinct id from systems where unique_id not in (SELECT distinct ase_id from reaction_system));".format(table))

        con.commit()
        con.close()
