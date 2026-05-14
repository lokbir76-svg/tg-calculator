import streamlit as st
import pandas as pd
import itertools
import io
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from itertools import combinations

# =============================================================================
# === КОНСТАНТЫ И СЛОВАРИ ===
# =============================================================================
MW_FA = {
    '14:0': 228.37, '16:0': 256.42, '18:0': 284.48, 
    '20:0': 312.53, '22:0': 340.58, '24:0': 368.64,
    '16:1': 254.41, '18:1': 282.46, '20:1': 310.51,
    '18:2': 280.45, '18:3': 278.43, '20:4': 304.47,
    '22:1': 336.55, '24:1': 364.60
}
GLYCEROL_BACKBONE_MW = 38.05

# =============================================================================
# === ФУНКЦИИ РАСЧЕТА ===
# =============================================================================

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
    """Прямой расчет: ЖК (масс.%) -> ТГ (масс.%)"""
    df_fa = df_input.copy()
    df_fa.columns = ['FA_Name', 'Mass_Percent']
    df_fa.dropna(inplace=True)
    df_fa['Mass_Percent'] = pd.to_numeric(df_fa['Mass_Percent'], errors='coerce')
    df_fa.dropna(subset=['Mass_Percent'], inplace=True)
    
    total_fa_mass = df_fa['Mass_Percent'].sum()
    if total_fa_mass > 0:
        df_fa['Mass_Percent'] = (df_fa['Mass_Percent'] / total_fa_mass) * 100
        
    df_fa['MW'] = df_fa['FA_Name'].apply(get_fa_mw)
    df_fa['Moles'] = df_fa['Mass_Percent'] / df_fa['MW']
    total_moles_fa = df_fa['Moles'].sum()
    df_fa['Mol_Fraction'] = df_fa['Moles'] / total_moles_fa if total_moles_fa > 0 else 0.0
    
    fa_dict = dict(zip(df_fa['FA_Name'].astype(str), df_fa['Mol_Fraction']))
    fa_list = list(fa_dict.keys())
    
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
    
    df_tg['MW_TG'] = df_tg['Triglyceride'].apply(calculate_tg_mw)
    df_tg['Mass_Contribution'] = df_tg['Mol_Fraction'] * df_tg['MW_TG']
    total_mass_contrib = df_tg['Mass_Contribution'].sum()
    
    df_tg['Mass_Percent'] = (df_tg['Mass_Contribution'] / total_mass_contrib * 100) if total_mass_contrib > 0 else 0.0
    df_tg['Mass_Percent'] = df_tg['Mass_Percent'].round(4)
    
    final_df = df_tg[['Triglyceride', 'Mass_Percent']].copy()
    final_df.sort_values(by='Mass_Percent', ascending=False, inplace=True)
    
    return final_df, df_fa

def calculate_fa_from_tg(df_tg_input):
    """Обратный расчет: ТГ (масс.%) -> ЖК (масс.%). Математически инвертирует process_data."""
    # Строго выбираем нужные столбцы
    if 'Triglyceride' not in df_tg_input.columns or 'Mass_Percent' not in df_tg_input.columns:
        # Фолбэк: берем первые две колонки, если заголовки отличаются
        df = df_tg_input.iloc[:, :2].copy()
        df.columns = ['Triglyceride', 'Mass_Percent']
    else:
        df = df_tg_input[['Triglyceride', 'Mass_Percent']].copy()
        
    df['Mass_Percent'] = pd.to_numeric(df['Mass_Percent'], errors='coerce')
    df = df.dropna(subset=['Mass_Percent'])
    if df.empty:
        return pd.DataFrame(columns=['FA_Name', 'Mass_Percent', 'Mol_Percent'])
        
    # Нормировка ТГ до 100%
    total_tg = df['Mass_Percent'].sum()
    df['Mass_Percent'] = (df['Mass_Percent'] / total_tg) * 100 if total_tg > 0 else 0.0
    
    # ТГ масс.% -> ТГ моль%
    df['MW_TG'] = df['Triglyceride'].apply(calculate_tg_mw)
    df['Moles_TG'] = df['Mass_Percent'] / df['MW_TG']
    total_moles_tg = df['Moles_TG'].sum()
    df['MolFrac_TG'] = df['Moles_TG'] / total_moles_tg if total_moles_tg > 0 else 0.0
    
    # Декомпозиция ТГ на ЖК (сбор мольных долей)
    fa_moles = {}
    for _, row in df.iterrows():
        tg_name = str(row['Triglyceride'])
        fas = tg_name.split('/')
        mol_frac_tg = row['MolFrac_TG']
        
        # Считаем сколько раз каждая кислота встречается в этом ТГ
        for fa in set(fas):
            fa_clean = fa.strip()
            count = fas.count(fa)
            # Вклад = мольная доля ТГ * (доля кислоты в молекуле ТГ)
            fa_moles[fa_clean] = fa_moles.get(fa_clean, 0.0) + mol_frac_tg * (count / 3.0)
            
    if not fa_moles:
        return pd.DataFrame(columns=['FA_Name', 'Mass_Percent', 'Mol_Percent'])
        
    df_fa = pd.DataFrame(list(fa_moles.items()), columns=['FA_Name', 'Mol_Percent'])
    
    # Нормировка ЖК мольных долей до 100%
    total_fa_mol = df_fa['Mol_Percent'].sum()
    df_fa['Mol_Percent'] = (df_fa['Mol_Percent'] / total_fa_mol * 100) if total_fa_mol > 0 else 0.0
    
    # ЖК моль% -> ЖК масс.%
    df_fa['MW_FA'] = df_fa['FA_Name'].apply(get_fa_mw)
    df_fa['Mass_Contrib'] = df_fa['Mol_Percent'] * df_fa['MW_FA']
    total_mass_fa = df_fa['Mass_Contrib'].sum()
    df_fa['Mass_Percent'] = (df_fa['Mass_Contrib'] / total_mass_fa * 100) if total_mass_fa > 0 else 0.0
    df_fa['Mass_Percent'] = df_fa['Mass_Percent'].round(3)
    
    return df_fa[['FA_Name', 'Mass_Percent', 'Mol_Percent']].sort_values('Mass_Percent', ascending=False)

# =============================================================================
# === ИНТЕРФЕЙС STREAMLIT ===
# =============================================================================

st.set_page_config(page_title="Калькулятор ТГ", layout="wide")
st.title("🧪 Расчет состава триглицеридов")
st.markdown("""
Загрузите Excel-файл с жирнокислотным составом (две колонки: Название ЖК и Массовая доля %).
Программа автоматически переведет данные в молярные доли, рассчитает распределение ТГ 
и вернет результат в массовых процентах.
""")

tab1, tab2, tab3, tab4 = st.tabs([
    "📊 Калькулятор", 
    "📝 Сравнение с экспериментом", 
    "💬 Обратная связь",
    "🔄 ТГ → ЖК (обратный)"
])

# --- Вкладка 1: Калькулятор ---
with tab1:
    uploaded_file = st.file_uploader("Выберите файл .xlsx или .csv", type=["xlsx", "csv"])
    
    if uploaded_file is not None:
        try:
            if uploaded_file.name.endswith('.csv'):
                df_input = pd.read_csv(uploaded_file, header=None)
            else:
                df_input = pd.read_excel(uploaded_file, header=None)
                
            if df_input.shape[1] < 2:
                st.error("Ошибка: В файле должно быть минимум 2 колонки.")
            else:
                df_input = df_input.iloc[:, :2]
                
                if st.button("🚀 Рассчитать"):
                    with st.spinner('Выполняется расчет...'):
                        result_df, fa_df = process_data(df_input)
                        st.success("Расчет успешно завершен!")
                        
                        t1, t2 = st.tabs(["📊 Результат (ТГ)", " Исходные данные (ЖК)"])
                        
                        with t1:
                            st.subheader("Массовое распределение триглицеридов")
                            st.dataframe(result_df, use_container_width=True)
                            
                            # === Расчет мольных процентов для экспорта ===
                            temp_df = result_df.copy()
                            temp_df['MW_TG'] = temp_df['Triglyceride'].apply(calculate_tg_mw)
                            temp_df['Moles_TG'] = temp_df['Mass_Percent'] / temp_df['MW_TG']
                            total_moles = temp_df['Moles_TG'].sum()
                            temp_df['Mol_Percent'] = (temp_df['Moles_TG'] / total_moles * 100).round(4) if total_moles > 0 else 0.0
                            
                            # Строгий порядок столбцов для совместимости с вкладкой 4
                            export_df = temp_df[['Triglyceride', 'Mass_Percent', 'Mol_Percent']].copy()
                            
                            st.subheader("📥 Скачать результаты")
                            col_dl1, col_dl2 = st.columns(2)
                            
                            with col_dl1:
                                st.download_button(
                                    label="📄 Скачать CSV",
                                    data=result_df.to_csv(index=False, encoding='utf-8-sig').encode('utf-8'),
                                    file_name='tg_distribution.csv',
                                    mime='text/csv'
                                )
                            
                            with col_dl2:
                                buffer = io.BytesIO()
                                with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                                    export_df.to_excel(writer, index=False, sheet_name='TG_Composition')
                                    fa_df.to_excel(writer, index=False, sheet_name='FA_Input')
                                buffer.seek(0)
                                st.download_button(
                                    label=" Скачать Excel (.xlsx)",
                                    data=buffer,
                                    file_name='tg_distribution_full.xlsx',
                                    mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                    help="Содержит столбцы: Triglyceride, Mass_Percent, Mol_Percent. Полностью совместим с вкладкой 'ТГ → ЖК'."
                                )
                            
                            st.caption("💡 Скачанный Excel-файл можно сразу загрузить во вкладку **'🔄 ТГ → ЖК'** для проверки обратного пересчета.")
                            
                            top_10 = result_df.head(10)
                            st.bar_chart(top_10.set_index('Triglyceride')['Mass_Percent'])

                        with t2:
                            st.subheader("Загруженный жирнокислотный состав")
                            st.dataframe(fa_df, use_container_width=True)
                            
        except Exception as e:
            st.error(f"Произошла ошибка: {e}")
            st.info("Убедитесь, что файл содержит два столбца: название кислоты и процент.")
    else:
        st.info("👆 Пожалуйста, загрузите файл выше, чтобы начать работу.")

# --- Вкладка 2: Сравнение с экспериментом ---
# --- Вкладка 2: Сравнение с экспериментом ---
with tab2:
    st.subheader("🔬 Сравнение структуры ТГ и состава ЖК")
    st.markdown("""
    Загрузите два файла:
    1. **Жирнокислотный состав (ЖКС)**: для расчета теоретического распределения ТГ
    2. **Экспериментальный профиль ТГ**: реальные данные (из масс-спектрометрии или хроматографии)
    """)
    
    col1, col2 = st.columns(2)
    with col1:
        fa_file_exp = st.file_uploader("1. Файл ЖКС (.xlsx/.csv)", type=["xlsx", "csv"], key="fa_comp")
    with col2:
        exp_file = st.file_uploader("2. Файл эксп. ТГ (.xlsx/.csv)", type=["xlsx", "csv"], key="exp_tg_comp")

    if exp_file and fa_file_exp:
        try:
            # Чтение ЖКС
            if fa_file_exp.name.endswith('.csv'):
                df_fa_input = pd.read_csv(fa_file_exp, header=None)
            else:
                df_fa_input = pd.read_excel(fa_file_exp, header=None)
            df_fa_input = df_fa_input.iloc[:, :2]
            
            # Расчет теоретического распределения ТГ
            result_theor, fa_processed = process_data(df_fa_input)
            
            # Чтение экспериментальных ТГ
            if exp_file.name.endswith('.csv'):
                df_exp_raw = pd.read_csv(exp_file)
            else:
                df_exp_raw = pd.read_excel(exp_file)
            
            # Обработка экспериментальных данных
            if df_exp_raw.shape[1] >= 2:
                # Проверяем, есть ли заголовки
                if 'Triglyceride' in df_exp_raw.columns or 'triglyceride' in str(df_exp_raw.columns[0]).lower():
                    df_exp_raw.columns = ['Triglyceride', 'Mass_Percent'] + [f'col_{i}' for i in range(2, df_exp_raw.shape[1])]
                else:
                    df_exp_raw.columns = ['Triglyceride', 'Mass_Percent'] + [f'col_{i}' for i in range(2, df_exp_raw.shape[1])]
                    
                df_exp_raw['Mass_Percent'] = pd.to_numeric(df_exp_raw['Mass_Percent'], errors='coerce')
                df_exp_raw = df_exp_raw.dropna(subset=['Mass_Percent', 'Triglyceride'])
            else:
                st.error("Файл ТГ должен иметь минимум 2 колонки.")
                st.stop()
            
            # === ДИАГНОСТИКА ===
            st.info(f"📊 **Теоретических ТГ:** {len(result_theor)} | **Экспериментальных ТГ:** {len(df_exp_raw)}")
            
            # Показываем топ-5 ТГ из каждого файла
            col_diag1, col_diag2 = st.columns(2)
            with col_diag1:
                st.caption("Топ-5 теоретических ТГ:")
                st.dataframe(result_theor.head(5), use_container_width=True)
            with col_diag2:
                st.caption("Топ-5 экспериментальных ТГ:")
                st.dataframe(df_exp_raw.head(5), use_container_width=True)
            
            # Находим общие ТГ
            common_tgs = set(result_theor['Triglyceride']) & set(df_exp_raw['Triglyceride'])
            st.caption(f"✅ **Общих ТГ:** {len(common_tgs)}")
            
            if len(common_tgs) < 2:
                st.warning("⚠️ Слишком мало общих триглицеридов для сравнения!")
                st.info("💡 **Возможные причины:**\n- Формат названий ТГ отличается (например, '16:0/18:1/18:2' vs '18:1/16:0/18:2')\n- В экспериментальном файле другие кислоты\n- Загрузите файл ТГ, полученный из вкладки 'Калькулятор'")
            else:
                # Фильтрация и выравнивание
                exp_filtered = df_exp_raw[df_exp_raw['Triglyceride'].isin(common_tgs)].sort_values('Mass_Percent', ascending=False).head(30)
                theor_filtered = result_theor[result_theor['Triglyceride'].isin(exp_filtered['Triglyceride'])]
                theor_aligned = theor_filtered.set_index('Triglyceride').reindex(exp_filtered['Triglyceride']).reset_index()
                
                # Функции для расчета вклада кислот
                def decompose_tg(tg_name):
                    if not isinstance(tg_name, str): return {}
                    parts = str(tg_name).split('/')
                    counts = {}
                    for p in parts:
                        counts[p.strip()] = counts.get(p.strip(), 0) + 1
                    return counts

                def calculate_acid_contribution(df_tg, percent_col='Mass_Percent'):
                    all_fas = sorted(set([fa for tg in df_tg['Triglyceride'] for fa in decompose_tg(tg).keys()]))
                    matrix = pd.DataFrame(index=df_tg['Triglyceride'].values, columns=all_fas, data=0.0)
                    for _, row in df_tg.iterrows():
                        tg_mass_pct = row[percent_col]
                        for fa, count in decompose_tg(row['Triglyceride']).items():
                            if fa in matrix.columns:
                                matrix.at[row['Triglyceride'], fa] = tg_mass_pct * (count / 3.0)
                    return matrix

                # Расчет матриц
                mat_exp = calculate_acid_contribution(exp_filtered, 'Mass_Percent')
                mat_theor = calculate_acid_contribution(theor_aligned, 'Mass_Percent')
                mat_diff = mat_exp.subtract(mat_theor, fill_value=0)
                
                # Визуализация
                st.subheader("🗺️ Карта отклонений (Эксперимент - Теория)")
                mat_diff_clean = mat_diff.fillna(0).astype(float)
                max_val = abs(mat_diff_clean).max().max()
                
                if max_val == 0 or np.isnan(max_val):
                    st.warning("⚠️ Все отклонения равны 0. Возможно, данные идентичны.")
                    max_val = 1
                
                fig_heat = px.imshow(
                    mat_diff_clean, 
                    color_continuous_scale='RdBu_r', 
                    range_color=[-max_val, max_val], 
                    aspect='auto',
                    labels=dict(x="Жирная кислота", y="Триглицерид", color="Отклонение")
                )
                fig_heat.update_layout(height=600, width=800)
                st.plotly_chart(fig_heat, use_container_width=True)
                
                # Статистика
                st.caption(f"Диапазон отклонений: [{mat_diff_clean.min().min():.4f}, {mat_diff_clean.max().max():.4f}]")
                
        except Exception as e:
            st.error(f"❌ Ошибка анализа: {e}")
            import traceback
            st.code(traceback.format_exc())
    else:
        st.info("👆 Загрузите оба файла для сравнения.")

# --- Вкладка 3: Обратная связь ---
with tab3:
    st.subheader("💬 Отзывы и предложения")
    st.write("Напишите нам!")
    with st.form("feedback_form"):
        feedback_text = st.text_area("Ваш отзыв:", height=150)
        submitted = st.form_submit_button("📤 Отправить")
    if submitted:
        if not feedback_text.strip():
            st.error("Введите текст!")
        else:
            st.success("✅ Спасибо! Отзыв принят.")

# --- Вкладка 4: Обратный пересчет ТГ → ЖК ---
with tab4:
    st.subheader("🔄 Обратный пересчет: ТГ → ЖК")
    st.markdown("""
    Загрузите файл с составом триглицеридов. Программа восстановит исходный жирнокислотный состав.
    Если файл получен из вкладки 'Калькулятор', результат должен совпасть с исходными данными (±0.01%).
    """)
    
    tg_file_reverse = st.file_uploader("Загрузите файл ТГ (.xlsx/.csv)", type=["xlsx", "csv"], key="tg_reverse")
    
    if tg_file_reverse is not None:
        try:
            df_tg_raw = pd.read_excel(tg_file_reverse) if not tg_file_reverse.name.endswith('.csv') else pd.read_csv(tg_file_reverse)
            
            # Защита от несовпадения столбцов
            if 'Triglyceride' not in df_tg_raw.columns or 'Mass_Percent' not in df_tg_raw.columns:
                df_tg_clean = df_tg_raw.iloc[:, :2].copy()
                df_tg_clean.columns = ['Triglyceride', 'Mass_Percent']
            else:
                df_tg_clean = df_tg_raw[['Triglyceride', 'Mass_Percent']].copy()
                
            if st.button("🔄 Выполнить обратный пересчет", key="btn_reverse"):
                with st.spinner('Пересчет...'):
                    result_fa = calculate_fa_from_tg(df_tg_clean)
                    if result_fa.empty:
                        st.error("Не удалось выделить кислоты. Проверьте формат названий ТГ (например, '16:0/18:1/18:2').")
                    else:
                        st.success("✅ Пересчет завершен!")
                        st.subheader("Восстановленный жирнокислотный состав")
                        st.dataframe(result_fa, use_container_width=True)
                        
                        # === КНОПКИ СКАЧИВАНИЯ (ДОБАВЛЕНО) ===
                        st.subheader("📥 Скачать результаты")
                        col_dl1, col_dl2 = st.columns(2)
                        
                        with col_dl1:
                            # Скачать как CSV
                            csv_fa = result_fa.to_csv(index=False, encoding='utf-8-sig').encode('utf-8')
                            st.download_button(
                                label="📄 Скачать как CSV",
                                data=csv_fa,
                                file_name='fa_composition_restored.csv',
                                mime='text/csv',
                            )
                        
                        with col_dl2:
                            # Скачать как Excel
                            buffer = io.BytesIO()
                            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                                result_fa.to_excel(writer, index=False, sheet_name='FA_Composition')
                            buffer.seek(0)
                            st.download_button(
                                label="📊 Скачать как Excel (.xlsx)",
                                data=buffer,
                                file_name='fa_composition_restored.xlsx',
                                mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                            )
                        # === КОНЕЦ КНОПОК ===
                        
                        top_10_fa = result_fa.head(10)
                        fig_pie = px.pie(top_10_fa, values='Mass_Percent', names='FA_Name', title='Топ-10 ЖК (масс.%)', hole=0.3)
                        st.plotly_chart(fig_pie, use_container_width=True)
                        
                        st.caption(f"Уникальных кислот: **{len(result_fa)}** | Сумма масс. %: **{result_fa['Mass_Percent'].sum():.2f}%**")
                        
        except Exception as e:
            st.error(f"Ошибка пересчета: {e}")
    else:
        st.info("👆 Загрузите файл ТГ для обратного пересчета.")
