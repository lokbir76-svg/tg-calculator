import streamlit as st
import pandas as pd
import itertools
import io # Для работы с файлами в памяти

# --- КОНСТАНТЫ И СЛОВАРИ (из вашего скрипта) ---
MW_FA = {
    '14:0': 228.37, '16:0': 256.42, '18:0': 284.48, 
    '20:0': 312.53, '22:0': 340.58, '24:0': 368.64,
    '16:1': 254.41, '18:1': 282.46, '20:1': 310.51,
    '18:2': 280.45, '18:3': 278.43, '20:4': 304.47,
    '22:1': 336.55, '24:1': 364.60
}
GLYCEROL_BACKBONE_MW = 38.05

# --- ФУНКЦИИ РАСЧЕТА (адаптированные из calc_tg.py) ---

def get_fa_mw(fa_name):
    if pd.isna(fa_name): return 0.0
    clean_name = str(fa_name).split()[0]
    return MW_FA.get(clean_name, 280.0)

def calculate_tg_mw(tg_name):
    if not isinstance(tg_name, str): return 0.0
    try:
        fas = tg_name.split('/')
        total_mw = GLYCEROL_BACKBONE_MW
        for fa in fas:
            total_mw += get_fa_mw(fa)
        return total_mw
    except Exception:
        return 0.0

def process_data(df_input):
    """Основная логика: Масс % ЖК -> Мол % ЖК -> Распределение ТГ -> Мол % ТГ -> Масс % ТГ"""
    
    # 1. Очистка и пересчет ЖК в молярные доли
    df_fa = df_input.copy()
    df_fa.columns = ['FA_Name', 'Mass_Percent'] # Предполагаем 2 колонки без заголовков или с ними
    df_fa.dropna(inplace=True)
    df_fa['Mass_Percent'] = pd.to_numeric(df_fa['Mass_Percent'], errors='coerce')
    df_fa.dropna(subset=['Mass_Percent'], inplace=True)
    
    # Нормировка суммы ЖК до 100%
    total_fa_mass = df_fa['Mass_Percent'].sum()
    if total_fa_mass > 0:
        df_fa['Mass_Percent'] = (df_fa['Mass_Percent'] / total_fa_mass) * 100
        
    df_fa['MW'] = df_fa['FA_Name'].apply(get_fa_mw)
    df_fa['Moles'] = df_fa['Mass_Percent'] / df_fa['MW']
    total_moles_fa = df_fa['Moles'].sum()
    df_fa['Mol_Fraction'] = df_fa['Moles'] / total_moles_fa
    
    fa_dict = dict(zip(df_fa['FA_Name'].astype(str), df_fa['Mol_Fraction']))
    fa_list = list(fa_dict.keys())
    
    # 2. Генерация ТГ и расчет их мольных долей
    tg_results = []
    for combo in itertools.combinations_with_replacement(fa_list, 3):
        unique_fas = set(combo)
        prob = 0.0
        if len(unique_fas) == 3:
            prob = 6 * fa_dict[combo[0]] * fa_dict[combo[1]] * fa_dict[combo[2]]
        elif len(unique_fas) == 2:
            counts = {x: combo.count(x) for x in unique_fas}
            double_fa = [k for k, v in counts.items() if v == 2][0]
            single_fa = [k for k, v in counts.items() if v == 1][0]
            prob = 3 * (fa_dict[double_fa]**2) * fa_dict[single_fa]
        else:
            prob = 1 * (fa_dict[combo[0]]**3)
            
        sorted_combo = sorted(combo)
        tg_name = "/".join(sorted_combo)
        tg_results.append({'Triglyceride': tg_name, 'Mol_Fraction': prob})
        
    df_tg = pd.DataFrame(tg_results)
    df_tg.sort_values(by='Mol_Fraction', ascending=False, inplace=True)
    
    # 3. Пересчет ТГ из мольных в массовые проценты
    df_tg['MW_TG'] = df_tg['Triglyceride'].apply(calculate_tg_mw)
    df_tg['Mass_Contribution'] = df_tg['Mol_Fraction'] * df_tg['MW_TG']
    total_mass_contrib = df_tg['Mass_Contribution'].sum()
    
    if total_mass_contrib > 0:
        df_tg['Mass_Percent'] = (df_tg['Mass_Contribution'] / total_mass_contrib) * 100
    else:
        df_tg['Mass_Percent'] = 0
        
    df_tg['Mass_Percent'] = df_tg['Mass_Percent'].round(3)
    
    # Итоговая таблица для вывода
    final_df = df_tg[['Triglyceride', 'Mass_Percent']].copy()
    final_df.sort_values(by='Mass_Percent', ascending=False, inplace=True)
    
    return final_df, df_fa # Возвращаем и результат ТГ, и обработанные данные ЖК

# --- ИНТЕРФЕЙС STREAMLIT ---

st.set_page_config(page_title="Калькулятор ТГ", layout="wide")
st.title("🧪 Расчет состава триглицеридов")
st.markdown("""
Загрузите Excel-файл с жирнокислотным составом (две колонки: Название ЖК и Массовая доля %).
Программа автоматически переведет данные в молярные доли, рассчитает распределение ТГ 
и вернет результат в массовых процентах.
""")

# Загрузка файла
uploaded_file = st.file_uploader("Выберите файл .xlsx или .csv", type=["xlsx", "csv"])

if uploaded_file is not None:
    try:
        # Чтение файла
        if uploaded_file.name.endswith('.csv'):
            df_input = pd.read_csv(uploaded_file, header=None)
        else:
            df_input = pd.read_excel(uploaded_file, header=None)
            
        # Проверка, что есть хотя бы 2 колонки
        if df_input.shape[1] < 2:
            st.error("Ошибка: В файле должно быть минимум 2 колонки (Название ЖК и Процент).")
        else:
            # Берем первые две колонки
            df_input = df_input.iloc[:, :2]
            
            if st.button("🚀 Рассчитать"):
                with st.spinner('Выполняется расчет...'):
                    result_df, fa_df = process_data(df_input)
                    
                    st.success("Расчет успешно завершен!")
                    
                    # Вкладки для отображения данных
                    tab1, tab2 = st.tabs(["📊 Результат (ТГ)", "📝 Исходные данные (ЖК)"])
                    
                    with tab1:
                        st.subheader("Массовое распределение триглицеридов")
                        st.dataframe(result_df, use_container_width=True)
                        
                        # Кнопка скачивания
                        csv = result_df.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label="📥 Скачать результат (CSV)",
                            data=csv,
                            file_name='tg_distribution_mass.csv',
                            mime='text/csv',
                        )
                        
                        # График топ-10 ТГ
                        top_10 = result_df.head(10)
                        st.bar_chart(top_10.set_index('Triglyceride')['Mass_Percent'])

                    with tab2:
                        st.subheader("Загруженный жирнокислотный состав")
                        st.dataframe(fa_df, use_container_width=True)
                        
    except Exception as e:
        st.error(f"Произошла ошибка при обработке файла: {e}")
        st.info("Убедитесь, что файл содержит два столбца: название кислоты (например, 18:2) и её процентное содержание.")

else:
    st.info("👆 Пожалуйста, загрузите файл выше, чтобы начать работу.")