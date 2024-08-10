import streamlit as st
from streamlit_ketcher import st_ketcher
import sys 

from ligify import __version__ as ligify_version

from ligify.format_results import format_results
from ligify.fetch_data import fetch_data
from ligify.predict.database.chemical_db import blast_chemical
from ligify.predict.pubchem import get_smiles, get_inchikey, get_name

# This is essentially the frontend component of the Ligify Web applications.

def run_ligify(chem, results, progress, chemical, filters):

    m_spacer1, metrics_col, m_spacer2 = results.container().columns((1,3,1))
    regulator_column, data_column = results.columns([1,3])

    if st.session_state.SUBMITTED:

        if chemical["smiles"] == None:
            data_column.subheader("Chemical input was not recognized. Please try a different input method.")

        else:

            # SMILES = str(chemical_smiles)
            # chem.image(f'http://hulab.rxnfinder.org/smi2img/{SMILES}/', width=200)

            regulators, metrics = fetch_data(chemical["InChiKey"], filters)

            select_spacerL, please_select, select_spacerR  = data_column.container().columns([1,2,1])     

            format_results(data_column, chemical["name"])

            regulator_column.header('')
            regulator_column.subheader('Sensor candidates')
            regulator_column.divider()

            # If no regulators are returned, suggest alternative queries
            if regulators == None:
                similar_chemicals = blast_chemical(chemical["smiles"], filters["max_alt_chems"])
                regulator_column.write("No regulators found")
                please_select.subheader("No associated reactions   :pensive:") 
                please_select.write("Consider an alternative query")   
                data_column.dataframe(similar_chemicals)
                
            # If regulators are returned, format display
            else:

                # Metrics data
                metrics_col.subheader("Search metrics")
                m_rhea, m_genes, m_filtered, m_operons, m_regs = metrics_col.columns(5)
                m_rhea.metric("Rhea reactions", metrics["RHEA Reactions"])
                m_genes.metric("Bacterial genes", metrics["Total genes"])
                m_filtered.metric("Filtered genes", metrics["Filtered genes"])
                m_operons.metric("Operons", metrics["Total operons"])
                m_regs.metric("Regulators", metrics["Total regulators"])
                metrics_col.divider()

                if not st.session_state.data:
                    please_select.subheader("Please select a regulator") 

                reg_acc_col, rank_col = regulator_column.columns((2,1))
                reg_acc_col.markdown("<h5>Regulator</h5>", unsafe_allow_html=True)
                rank_col.markdown("<h5>Rank</h5>", unsafe_allow_html=True)

                for i in range(0, len(regulators)):
                    name = "var"+str(i)
                    rank = regulators[i]["rank"]["rank"]
                    color = regulators[i]["rank"]["color"]
                    rank_col.markdown(f"<p style='font-size:20px; font-weight: 600; color: {color};'>{rank}</p>", unsafe_allow_html=True)
                    name = reg_acc_col.form_submit_button(regulators[i]['refseq'])
                    if name:
                        st.session_state.data = regulators[i]
                        st.experimental_rerun()
