"""
標本資料庫品質稽核 (Specimen Database Quality Audit)
====================================================
針對 NCHU_Pringles 蜘蛛標本資料庫的規則式資料品質檢核。

設計理念
--------
建立一個「規則註冊表」框架
每條規則是一個獨立函式，回傳違規紀錄與嚴重度。
新增檢核項目只需要註冊新函式，不需要改動主流程。

對應實際運用的資料品質框架（如 Great Expectations），
也對應到臨床試驗中的 edit check / data validation plan。

Author: 董亭妤 (Ting-Yu Tung)
"""

import pandas as pd

# ---------------------------------------------------------------------------
# 設定：領域規則常數
# ---------------------------------------------------------------------------

# 台灣本島與離島的合理座標範圍
TAIWAN_BOUNDS = {"lat": (21.5, 25.5), "lon": (119.5, 122.5)}

# 儲位編碼的設計慣例：每個儲位盒對應一個屬
STORAGE_GENUS_RULE = {
    "Pringles_01": "Aoaraneus",
    "Pringles_02": "Bijoaraneus",
    # Pringles_03 為混合盒，不設限
}

# 定序工作流程的先後順序：必須先萃取 DNA，才能擴增各基因片段
WORKFLOW_ORDER = ["DNA", "COI", "28S"]

# 未描述種的命名模式——這些沒有中文俗名是「結構性缺失」，不是資料錯誤
UNDESCRIBED_PATTERN = r"(?:sp\.?\d*$|cf\.)"

# 經去識別化處理、已一般化至縣市層級的地名。
# 這些紀錄的「一個地名對應多組座標」是刻意設計的結果，不是資料錯誤，
# 因此在 GEO-02 中必須排除，否則會與去識別化政策互相衝突。
GENERALISED_LOCALITY_PATTERN = r"(?:County|City), (?:Taiwan|India|Philippines|Japan)$"


# ---------------------------------------------------------------------------
# 規則註冊表
# ---------------------------------------------------------------------------

RULES = []


def rule(code, description, severity):
    """裝飾器：將檢核函式註冊進規則表。"""

    def wrapper(func):
        RULES.append(
            {"code": code, "description": description,
             "severity": severity, "func": func}
        )
        return func

    return wrapper


# --- 值域與格式 -------------------------------------------------------------

@rule("GEO-01", "台灣樣本座標落在合理範圍外", "HIGH")
def check_coordinates(df):
    tw = df[df["Country"] == "Taiwan"]
    lat_lo, lat_hi = TAIWAN_BOUNDS["lat"]
    lon_lo, lon_hi = TAIWAN_BOUNDS["lon"]
    mask = (
        tw["Latitude"].between(lat_lo, lat_hi) &
        tw["Longitude"].between(lon_lo, lon_hi)
    )
    return tw[~mask & tw["Latitude"].notna()][
        ["ABARA_code", "Locality", "Latitude", "Longitude"]
    ]

@rule("GEO-02", "同一地點名稱對應座標不一致", "MEDIUM")
def check_locality(df):
    """
    檢查同一 Locality 字串是否指向一致的座標。

    排除已一般化的地名：這些紀錄為保護未描述分類群而刻意降低了地點解析度，
    「一個縣市對應多組座標」是政策的預期結果，不是資料品質問題。
    若不排除，去識別化會使本規則產生大量偽陽性——資料保護措施與
    資料品質規則之間的這類衝突，需要在規則層級明確處理。
    """
    named = df.dropna(subset=["Latitude", "Longitude"])
    named = named[~named["Locality"].str.contains(
        GENERALISED_LOCALITY_PATTERN, regex=True, na=False)]
    counts = named.groupby("Locality")[["Latitude", "Longitude"]].nunique()
    conflicted = counts[(counts["Latitude"] > 1) | (counts["Longitude"] > 1)].index
    return named[named["Locality"].isin(conflicted)][
        ["ABARA_code", "Locality", "Latitude", "Longitude"]
    ].drop_duplicates(subset=["Locality", "Latitude", "Longitude"])


@rule("FMT-01", "採集日欄位遭試算表軟體自動轉換為日期格式", "HIGH")
def check_day_corruption(df):
    """
    原始資料日期欄中的『8-16』代表 8 到 16 日的採集期間，
    但 Excel/Sheets 會將其自動判讀為日期並改寫成 08/16/YYYY。
    這是 Ziemann et al. (2016) 指出的基因名稱被 Excel 破壞的同類問題。
    """
    day = df["Collection_Day"].astype(str)
    mask = day.str.contains("/", na=False)
    return df[mask][
        ["ABARA_code", "Collection_Year", "Collection_Month", "Collection_Day"]
    ]


# --- 類別編碼一致性 ---------------------------------------------------------

@rule("STD-01", "國別欄位的分類層級不一致", "MEDIUM")
def check_country_granularity(df):
    """沖繩 (Okinawa) 雖為離島但仍屬於日本的行政區，與 Japan 並列會造成分組統計錯誤。"""
    subregions = {"Okinawa"}
    return df[df["Country"].isin(subregions)][
        ["ABARA_code", "Country", "Locality"]
    ]


@rule("STD-02", "同一學名對應到多個中文俗名", "MEDIUM")
def check_vernacular_consistency(df):
    named = df.dropna(subset=["Chinese_name"])
    counts = named.groupby("Scientific_name")["Chinese_name"].nunique()
    conflicted = counts[counts > 1].index
    return named[named["Scientific_name"].isin(conflicted)][
        ["ABARA_code", "Scientific_name", "Chinese_name"]
    ].drop_duplicates(subset=["Scientific_name", "Chinese_name"])


@rule("STD-03", "以字串哨兵值表示缺失", "LOW")
def check_sentinel_values(df):
    """'Missing' 混入類別欄位，會被當成一個真實的儲位類別。"""
    sentinels = {"Missing", "N/A", "NA", "unknown", "-"}
    hits = []
    for col in df.select_dtypes(include=["object", "str"]).columns:
        m = df[col].isin(sentinels)
        if m.any():
            out = df[m][["ABARA_code", col]].copy()
            out.columns = ["ABARA_code", "value"]
            out.insert(1, "column", col)
            hits.append(out)
    return pd.concat(hits) if hits else pd.DataFrame()


# --- 業務規則與參照完整性 ---------------------------------------------------

@rule("REF-01", "標本儲位與所屬屬別不符", "MEDIUM")
def check_storage_rule(df):
    violations = []
    for storage, genus in STORAGE_GENUS_RULE.items():
        v = df[(df["Storage"] == storage) & (df["Genus"] != genus)]
        violations.append(v)
    out = pd.concat(violations)
    return out[["ABARA_code", "Genus", "Species", "Storage"]]


@rule("SEQ-01", "定序流程跳階：下游步驟已完成但上游未完成", "MEDIUM")
def check_workflow_order(df):
    """
    正常流程為 DNA 萃取 -> COI 擴增 -> 28S 擴增。
    若 28S 有紀錄而 COI 無，代表流程紀錄漏填或實際跳過了步驟。
    """
    rows = []
    for i in range(1, len(WORKFLOW_ORDER)):
        downstream = WORKFLOW_ORDER[i]
        upstream = WORKFLOW_ORDER[i - 1]
        m = df[downstream].notna() & df[upstream].isna()
        if m.any():
            out = df[m][["ABARA_code"] + WORKFLOW_ORDER].copy()
            out["missing_step"] = upstream
            out["completed_step"] = downstream
            rows.append(out)
    return pd.concat(rows) if rows else pd.DataFrame()


@rule("SEQ-02", "完成狀態符號在各欄位間不一致", "LOW")
def check_symbol_vocabulary(df):
    """
    DNA 與 28S 使用 '●'，COI 只使用 '◐'（代表非雙向定序）。
    同一個『已完成』的語意用了兩種符號，會讓程式化解析產生歧義。
    """
    vocab = {c: set(df[c].dropna().unique()) for c in WORKFLOW_ORDER}
    all_symbols = set().union(*vocab.values())
    inconsistent = any(v != all_symbols for v in vocab.values())
    if not inconsistent:
        return pd.DataFrame()
    return pd.DataFrame(
        [{"column": c, "symbols_used": ", ".join(sorted(s))}
         for c, s in vocab.items()]
    )


# --- 缺失值 -----------------------------------------------------------------

@rule("MIS-01", "隨機缺失的中文俗名（同種其他個體已有紀錄）", "LOW")
def check_random_missing_vernacular(df):
    """
    未描述種 (sp.1, sp.2, cf.) 沒有中文名屬於『結構性缺失』，是正常的。
    但若同一學名的其他標本已有中文名，單筆缺漏就是『隨機缺失』，應補齊。
    區分這兩者是缺失值處理的第一步。
    """
    described = df[~df["Scientific_name"].str.contains(
        UNDESCRIBED_PATTERN, regex=True, na=False)]
    has_name = described.dropna(subset=["Chinese_name"])
    known = set(has_name["Scientific_name"])
    m = described["Chinese_name"].isna() & described["Scientific_name"].isin(known)
    return described[m][["ABARA_code", "Scientific_name", "Chinese_name"]]


# ---------------------------------------------------------------------------
# 執行與報告
# ---------------------------------------------------------------------------

def missing_value_summary(df):
    """區分結構性缺失與需處理的缺失。"""
    total = len(df)
    rows = []
    for col in df.columns:
        n = int(df[col].isna().sum())
        if n:
            rows.append({"column": col, "n_missing": n,
                         "pct": round(100 * n / total, 1)})
    return pd.DataFrame(rows).sort_values("n_missing", ascending=False)


def run_audit(path):
    df = pd.read_csv(path)

    print("=" * 68)
    print("標本資料庫品質稽核報告")
    print("=" * 68)
    print(f"資料來源：{path}")
    print(f"資料規模：{df.shape[0]} 筆標本 × {df.shape[1]} 個欄位")
    print(f"分類涵蓋：{df['Genus'].nunique()} 屬 / "
          f"{df['Scientific_name'].nunique()} 個分類單元")
    print(f"地理涵蓋：{df['Country'].nunique()} 個國別標記")
    print(f"時間跨度：{df['Collection_Year'].min()}–{df['Collection_Year'].max()}")

    print("\n" + "-" * 68)
    print("一、缺失值概況")
    print("-" * 68)
    print(missing_value_summary(df).to_string(index=False))

    print("\n" + "-" * 68)
    print("二、規則檢核結果")
    print("-" * 68)

    summary = []
    for r in RULES:
        result = r["func"](df)
        n = len(result)
        summary.append({"code": r["code"], "severity": r["severity"],
                        "description": r["description"], "n_violations": n})
        if n:
            print(f"\n[{r['severity']}] {r['code']} — {r['description']}"
                  f"  ({n} 筆)")
            print(result.to_string(index=False))

    print("\n" + "-" * 68)
    print("三、稽核總表")
    print("-" * 68)
    s = pd.DataFrame(summary)
    print(s.to_string(index=False))
    print(f"\n合計 {int(s['n_violations'].sum())} 筆違規，"
          f"涉及 {int((s['n_violations'] > 0).sum())} 條規則。")

    return df, s


if __name__ == "__main__":
    run_audit("data/NCHU_Pringles.csv")
