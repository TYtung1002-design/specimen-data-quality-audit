"""
未描述分類群的地點資料去識別化 (Locality De-identification)
==========================================================
在公開標本資料庫前，降低未描述分類群（新種、新紀錄種）採集地點的解析度。

為什麼需要這個步驟
------------------
未正式發表的新種，其模式產地在論文出版前不宜公開；部分稀有類群也有
採集壓力的疑慮。但直接刪除整筆紀錄會讓資料集失去分析價值。

去識別化的目標不是「移除資訊」，而是「降低解析度到足以保護、
但仍可供分析」的程度——這與醫療資料的去識別化是同一個問題。

設計要點
--------
1. **同時處理直接識別欄位與準識別欄位。**
   只模糊座標而保留精確地名（例如「某某停車場」）是無效的，
   因為地名可反查出比模糊化後座標更高的精度。
   這對應到 HIPAA Safe Harbor 要求移除 18 類欄位而不只是姓名的理由。

2. **公開版是衍生資料，不是原始資料。**
   本腳本從完整資料產生公開版；完整資料不進入版本控制。

3. **無法自動處理的紀錄要進人工審查佇列，不能靜默放行。**

Author: 董亭妤 (Ting-Yu Tung)
"""

import re
import pandas as pd

# ---------------------------------------------------------------------------
# 去識別化政策設定
# ---------------------------------------------------------------------------

# 需要保護的分類群：尚未正式描述或未鑑定到種的類群
# 若後續某些類群已發表，可自此移除
SENSITIVE_TAXA_PATTERN = r"(?:\bsp\.?\s*\d*$|\bcf\.)"

# 座標保留小數位數。2 位約等於 1.1 公里的緯度解析度
COORD_DECIMALS = 2

# 海拔分組間距（公尺）。海拔在山區是很強的準識別欄位
ELEVATION_BIN = 100

# 地名一般化的目標層級：縣市
# 鍵為可在原始地名中比對到的樣式，值為一般化後的標準寫法
COUNTY_LOOKUP = [
    (r"屏東|Pingtung|港口溪|龍鑾潭|恆春",      "Pingtung County, Taiwan"),
    (r"嘉義|Chiayi",                "Chiayi County, Taiwan"),
    (r"南投|Nantou",                "Nantou County, Taiwan"),
    (r"台東|臺東|Taitung",           "Taitung County, Taiwan"),
    (r"宜蘭|Yilan|Ilan",            "Yilan County, Taiwan"),
    (r"花蓮|Hualien",               "Hualien County, Taiwan"),
    (r"苗栗|Miaoli",                "Miaoli County, Taiwan"),
    (r"雲林|Yunlin",                "Yunlin County, Taiwan"),
    (r"新竹|Hsinchu",               "Hsinchu County, Taiwan"),
    (r"新北|New Taipei",            "New Taipei City, Taiwan"),
    (r"陽明山|北投|大安|Taipei",      "Taipei City, Taiwan"),
    (r"台中|臺中|Taichung",          "Taichung City, Taiwan"),
    (r"高雄|Kaohsiung",             "Kaohsiung City, Taiwan"),
    (r"台南|臺南|Tainan",            "Tainan City, Taiwan"),
    (r"基隆|Kee ?Lung|Keelung",      "Keelung City, Taiwan"),
    (r"Tamil Nadu",                 "Tamil Nadu, India"),
    (r"Luzon|Apayao",               "Luzon, Philippines"),
    (r"Okinawa|沖繩",                "Okinawa, Japan"),
]

# 自動比對失敗時的暫代值——這些紀錄會被列入人工審查清單
UNRESOLVED_LABEL = "LOCALITY_WITHHELD_PENDING_REVIEW"


# ---------------------------------------------------------------------------
# 去識別化步驟
# ---------------------------------------------------------------------------

def flag_sensitive(df):
    """標記需要保護的紀錄。"""
    return df["Scientific_name"].str.contains(
        SENSITIVE_TAXA_PATTERN, regex=True, na=False)


def generalise_locality(text):
    """將精確地名一般化至縣市層級；無法判定者回傳暫代值。"""
    if pd.isna(text):
        return text
    for pattern, label in COUNTY_LOOKUP:
        if re.search(pattern, str(text)):
            return label
    return UNRESOLVED_LABEL


def deidentify(df):
    """回傳 (公開版資料, 處理摘要, 待人工審查清單)。"""
    out = df.copy()
    mask = flag_sensitive(out)

    # 1. 座標降低解析度
    for col in ["Latitude", "Longitude"]:
        out.loc[mask, col] = out.loc[mask, col].round(COORD_DECIMALS)

    # 2. 海拔分組
    out.loc[mask, "Elve"] = (
        (out.loc[mask, "Elve"] / ELEVATION_BIN).round() * ELEVATION_BIN
    )

    # 3. 地名一般化
    out.loc[mask, "Locality"] = out.loc[mask, "Locality"].apply(generalise_locality)

    pending = out[mask & (out["Locality"] == UNRESOLVED_LABEL)]

    summary = {
        "總筆數": len(df),
        "受保護紀錄": int(mask.sum()),
        "受保護分類群": sorted(df.loc[mask, "Scientific_name"].unique()),
        "座標已降解析度": int((mask & df["Latitude"].notna()).sum()),
        "地名已一般化": int((mask & df["Locality"].notna()).sum()),
        "待人工審查": len(pending),
    }
    return out, summary, pending


def report_precision_loss(before, after, mask):
    """量化去識別化造成的定位誤差，用於評估保護強度是否足夠。"""
    b = before.loc[mask & before["Latitude"].notna()]
    a = after.loc[mask & after["Latitude"].notna()]
    dlat = (a["Latitude"] - b["Latitude"]).abs() * 111.0
    dlon = (a["Longitude"] - b["Longitude"]).abs() * 101.0
    return {
        "最大南北位移_km": round(float(dlat.max()), 3),
        "最大東西位移_km": round(float(dlon.max()), 3),
        "平均位移_km": round(float(((dlat**2 + dlon**2) ** 0.5).mean()), 3),
    }


def main(src, dst):
    df = pd.read_csv(src)
    mask = flag_sensitive(df)
    out, summary, pending = deidentify(df)

    print("=" * 68)
    print("地點資料去識別化報告")
    print("=" * 68)
    for k, v in summary.items():
        if isinstance(v, list):
            print(f"{k}：")
            for item in v:
                print(f"  - {item}")
        else:
            print(f"{k}：{v}")

    print("\n定位精度損失（受保護紀錄）")
    for k, v in report_precision_loss(df, out, mask).items():
        print(f"  {k}：{v}")

    if len(pending):
        print("\n" + "!" * 68)
        print("以下紀錄無法自動判定縣市，需人工確認後補上：")
        print("!" * 68)
        print(pending[["ABARA_code", "Scientific_name"]].to_string(index=False))
        print("\n原始地名：")
        for code in pending["ABARA_code"]:
            print(f"  {code}: {df.loc[df.ABARA_code == code, 'Locality'].iloc[0]}")

    out.to_csv(dst, index=False)
    print(f"\n公開版已輸出：{dst}")
    return out


if __name__ == "__main__":
    main("data/NCHU_Pringles_full.csv", "data/NCHU_Pringles.csv")
